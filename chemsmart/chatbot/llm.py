"""LLM interface for chemsmart chatbot.

This module provides an LLM class that can work with:
1. OpenAI API (set OPENAI_API_KEY environment variable)
2. OpenAI-compatible APIs (set OPENAI_API_BASE and OPENAI_API_KEY)
3. Offline mode with template-based response generation
"""

import json
import logging
import os
import re

from . import COMMANDS_DATASET

logger = logging.getLogger(__name__)

# System prompt for the LLM
SYSTEM_PROMPT = """You are ChemSmart Assistant, an expert in computational chemistry \
and quantum chemistry simulations. Your task is to help users create and submit \
chemistry simulation jobs using the chemsmart command-line interface.

When a user describes a job they want to run, you should generate the appropriate \
chemsmart command.

Available job types for Gaussian:
- opt: Geometry optimization
- ts: Transition state search
- irc: Intrinsic reaction coordinate
- sp: Single-point energy calculation
- scan: Relaxed PES scan
- modred: Modredundant frozen coordinates optimization
- nci: Non-covalent interaction
- resp: RESP charges fitting
- crest: Optimize conformers from CREST output
- saopt: Optimize structures from MD trajectory
- dias: Distortion-interaction/activation-strain analysis
- userjob: User-defined custom job
- com: Run pre-prepared Gaussian input directly

Command format:
chemsmart sub -s <server> gaussian -p <project> -f <input_file> [options] <job_type>

Common options:
- -s/--server: Server name (e.g., shared, local)
- -p/--project: Project settings name (e.g., test)
- -f/--file: Input file path
- -c/--charge: System charge (default: from input file)
- -m/--multiplicity: System multiplicity (default: from input file)
- -l/--label: Custom output file name

When you generate a command, output ONLY the command starting with "chemsmart" \
without any explanation. If you cannot generate a command, explain why.
"""

# Pre-compiled regex patterns for performance
_FILE_PATTERN = re.compile(r"(\S+\.(com|xyz|log|gjf|inp|out|traj|db))")
_SERVER_PATTERN = re.compile(r"server\s+(\w+)")
_PROJECT_PATTERN = re.compile(r"project\s+(\w+)")
_FILE_KEYWORD_PATTERN = re.compile(r"file\s+(\S+)")

# Pre-compiled job type patterns
_JOB_TYPE_PATTERNS = [
    (re.compile(r"\birc\b"), "irc"),
    (re.compile(r"\bts\b"), "ts"),
    (re.compile(r"\bsp\b"), "sp"),
    (re.compile(r"\bscan\b"), "scan"),
    (re.compile(r"\bmodred\b"), "modred"),
    (re.compile(r"\bfrozen\b"), "modred"),
    (re.compile(r"\bnci\b"), "nci"),
    (re.compile(r"\bresp\b"), "resp"),
    (re.compile(r"\bcrest\b"), "crest"),
    (re.compile(r"\bconformer"), "crest"),
    (re.compile(r"\bdias\b"), "dias"),
    (re.compile(r"\bdistortion"), "dias"),
    (re.compile(r"\boptimiz"), "opt"),
    (re.compile(r"\bopt\b"), "opt"),
]


class LLM:
    """LLM interface that supports OpenAI API or offline mode."""

    def __init__(self, model=None, api_key=None, api_base=None):
        """Initialize the LLM.

        Args:
            model: Model name (default: gpt-4o-mini or from OPENAI_MODEL env)
            api_key: OpenAI API key (default: from OPENAI_API_KEY env)
            api_base: OpenAI API base URL (default: from OPENAI_API_BASE env)
        """
        self.model = model or os.environ.get("OPENAI_MODEL", "gpt-4o-mini")
        self.api_key = api_key or os.environ.get("OPENAI_API_KEY")
        self.api_base = api_base or os.environ.get(
            "OPENAI_API_BASE", "https://api.openai.com/v1"
        )
        self.client = None
        self.offline_mode = False

        # Load training data for offline mode
        self._load_training_data()

        # Try to initialize OpenAI client
        if self.api_key:
            try:
                from openai import OpenAI

                self.client = OpenAI(
                    api_key=self.api_key, base_url=self.api_base
                )
                logger.info(f"LLM initialized with model: {self.model}")
            except ImportError:
                logger.warning(
                    "openai package not installed. Running in offline mode."
                )
                self.offline_mode = True
        else:
            logger.warning(
                "OPENAI_API_KEY not set. Running in offline mode."
            )
            self.offline_mode = True

    def _load_training_data(self):
        """Load training data for offline pattern matching."""
        self.training_data = []
        try:
            with open(COMMANDS_DATASET, "r") as f:
                data = json.load(f)
                self.training_data = data.get("training_data", [])
        except (FileNotFoundError, json.JSONDecodeError) as e:
            logger.warning(f"Could not load training data: {e}")

    def _offline_response(self, user_input):
        """Generate response using pattern matching in offline mode."""
        user_lower = user_input.lower()

        # First, try to construct command from explicit keywords in user input
        # This is more reliable than pattern matching
        constructed_cmd = self._construct_command_from_keywords(user_input)

        # If we found explicit patterns (file, server, project), use constructed
        # command as it will be more accurate to user's intent
        has_explicit_params = (
            _FILE_PATTERN.search(user_lower)
            or _SERVER_PATTERN.search(user_lower)
            or _PROJECT_PATTERN.search(user_lower)
        )

        if has_explicit_params:
            return constructed_cmd

        # Try to match with training data for more complex examples
        best_match = None
        best_score = 0

        for entry in self.training_data:
            desc = entry.get("description", "").lower()
            # Simple keyword matching score
            score = sum(1 for word in user_lower.split() if word in desc)
            if score > best_score:
                best_score = score
                best_match = entry

        # Only use training data match if score is high enough
        if best_match and best_score >= 3:
            return best_match.get("command", "")

        # Fallback to constructed command
        return constructed_cmd

    def _construct_command_from_keywords(self, user_input):
        """Construct a command from keywords in user input."""
        user_lower = user_input.lower()

        # Detect job type - check more specific patterns first
        job_type = "opt"  # default

        # Check for multi-word patterns first (more specific)
        multi_word_mappings = [
            ("transition state", "ts"),
            ("single point", "sp"),
            ("single-point", "sp"),
            ("intrinsic reaction", "irc"),
            ("pes scan", "scan"),
            ("non-covalent", "nci"),
        ]

        for pattern, job in multi_word_mappings:
            if pattern in user_lower:
                job_type = job
                break
        else:
            # Check for single-word patterns with pre-compiled regexes
            for pattern, job in _JOB_TYPE_PATTERNS:
                if pattern.search(user_lower):
                    job_type = job
                    break

        # Detect server using pre-compiled pattern
        server = "shared"  # default
        server_match = _SERVER_PATTERN.search(user_lower)
        if server_match:
            server = server_match.group(1)

        # Detect project using pre-compiled pattern
        project = "test"  # default
        project_match = _PROJECT_PATTERN.search(user_lower)
        if project_match:
            project = project_match.group(1)

        # Detect input file using pre-compiled patterns
        input_file = "input.com"  # default
        file_match = _FILE_KEYWORD_PATTERN.search(user_lower)
        if file_match:
            input_file = file_match.group(1)
        else:
            # Look for common file extensions
            file_ext_match = _FILE_PATTERN.search(user_lower)
            if file_ext_match:
                input_file = file_ext_match.group(1)

        return (
            f"chemsmart sub -s {server} gaussian "
            f"-p {project} -f {input_file} {job_type}"
        )

    def generate_response(self, prompt):
        """Generate a response from the LLM.

        Args:
            prompt: User prompt

        Returns:
            Generated response string
        """
        if self.offline_mode or not self.client:
            return self._offline_response(prompt)

        try:
            response = self.client.chat.completions.create(
                model=self.model,
                messages=[
                    {"role": "system", "content": SYSTEM_PROMPT},
                    {"role": "user", "content": prompt},
                ],
                max_tokens=256,
                temperature=0.3,
            )
            return response.choices[0].message.content.strip()
        except Exception as e:
            logger.error(f"LLM API error: {e}")
            # Fallback to offline mode
            return self._offline_response(prompt)

    def chat(self, messages):
        """Have a multi-turn conversation with the LLM.

        Args:
            messages: List of message dicts with 'role' and 'content'

        Returns:
            Generated response string
        """
        if self.offline_mode or not self.client:
            # In offline mode, just use the last user message
            for msg in reversed(messages):
                if msg.get("role") == "user":
                    return self._offline_response(msg.get("content", ""))
            return "I couldn't understand your request. Please try again."

        try:
            full_messages = [{"role": "system", "content": SYSTEM_PROMPT}]
            full_messages.extend(messages)

            response = self.client.chat.completions.create(
                model=self.model,
                messages=full_messages,
                max_tokens=512,
                temperature=0.3,
            )
            return response.choices[0].message.content.strip()
        except Exception as e:
            logger.error(f"LLM API error: {e}")
            # Fallback to offline mode
            for msg in reversed(messages):
                if msg.get("role") == "user":
                    return self._offline_response(msg.get("content", ""))
            return "I couldn't understand your request. Please try again."
