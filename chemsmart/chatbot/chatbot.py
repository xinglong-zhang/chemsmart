"""ChemSmart Chatbot for interactive command generation and execution."""

import re

from .executor import CommandExecutor
from .llm import LLM
from .logger import log


class Chatbot:
    """Interactive chatbot for generating and executing chemsmart commands."""

    def __init__(self, auto_execute=False):
        """Initialize the chatbot.

        Args:
            auto_execute: If True, automatically execute generated commands
        """
        self.llm = LLM()
        self.executor = CommandExecutor()
        self.auto_execute = auto_execute
        self.conversation_history = []

    def _is_chemsmart_command(self, text):
        """Check if text contains a chemsmart command."""
        return "chemsmart" in text.lower()

    def _extract_command(self, response):
        """Extract chemsmart command from LLM response."""
        # Look for command on its own line or after EXECUTE:
        lines = response.split("\n")
        for line in lines:
            line = line.strip()
            if line.startswith("EXECUTE:"):
                return line.replace("EXECUTE:", "").strip()
            if line.startswith("chemsmart "):
                return line
            # Handle markdown code blocks
            if line.startswith("```"):
                continue
            # Look for command in backticks
            match = re.search(r"`(chemsmart[^`]+)`", line)
            if match:
                return match.group(1)
        return None

    def process_message(self, user_input):
        """Process a single user message and return response.

        Args:
            user_input: User's message

        Returns:
            Tuple of (response_text, command, command_output)
        """
        # Add to conversation history
        self.conversation_history.append(
            {"role": "user", "content": user_input}
        )

        # Get LLM response
        response = self.llm.chat(self.conversation_history)

        # Add assistant response to history
        self.conversation_history.append(
            {"role": "assistant", "content": response}
        )

        # Check if response contains a command
        command = self._extract_command(response)
        command_output = None

        if command and self.auto_execute:
            command_output = self.executor.execute(command)

        return response, command, command_output

    def chats(self):
        """Start an interactive chat session."""
        log("ChemSmart Chatbot started.")
        log("I can help you create and submit chemistry simulation jobs.")
        log("Describe what job you want to run, and I'll generate the command.")
        log("Type 'exit' or 'quit' to end the session.")
        log("Type 'execute' to run the last generated command.")
        log("Type 'clear' to clear conversation history.")
        log("-" * 50)

        last_command = None

        while True:
            try:
                user_input = input("\nYou: ").strip()
            except (EOFError, KeyboardInterrupt):
                log("\nChatbot session ended.")
                break

            if not user_input:
                continue

            if user_input.lower() in ("exit", "quit"):
                log("Chatbot session ended.")
                break

            if user_input.lower() == "clear":
                self.conversation_history = []
                log("Conversation history cleared.")
                continue

            if user_input.lower() == "execute":
                if last_command:
                    log(f"Executing: {last_command}")
                    output = self.executor.execute(last_command)
                    log(f"Output:\n{output}")
                else:
                    log("No command to execute. Generate a command first.")
                continue

            # Process the message
            response, command, command_output = self.process_message(user_input)

            print(f"\nAssistant: {response}")

            if command:
                last_command = command
                print(f"\nGenerated command: {command}")

                if command_output is not None:
                    print(f"\nCommand output:\n{command_output}")
                else:
                    print(
                        "\nType 'execute' to run this command, "
                        "or ask me to modify it."
                    )

    def generate_command(self, description):
        """Generate a command from a natural language description.

        Args:
            description: Natural language description of the job

        Returns:
            Generated chemsmart command string, or None if failed
        """
        response = self.llm.generate_response(description)
        return self._extract_command(response)

    def train(self):
        """Placeholder for training functionality.

        Note: Training requires additional dependencies
        (transformers, torch, datasets) and is not supported
        in the basic chatbot mode.
        """
        try:
            from .llm_trainer import ChemSmartTrainer

            trainer = ChemSmartTrainer()
            trainer.preprocess_data()
            trainer.train()
        except ImportError as e:
            log(f"Training requires additional dependencies: {e}")
            log("Install with: pip install transformers torch datasets")


if __name__ == "__main__":
    Chatbot().chats()
