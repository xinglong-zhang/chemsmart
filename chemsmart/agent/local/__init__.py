"""Local v13.1 model adapter for chemsmart agent.

Wraps the Smilesjs/chemsmart-qwen2.5-coder-3b-instruct-v13_1 model so that natural-language
queries can be served through the same ``chemsmart agent ask`` flow used for
OpenAI/Anthropic providers, without exposing remote API keys.

Submodules:
    loader        Load the merged v13.1 model on CUDA, MPS, or CPU.
    generator     Greedy decode wrapper that returns compact SPEC JSON.
    adapter       Convert compact SPEC to SynthesisSession schema.
    server        OpenAI-compatible /v1/chat/completions FastAPI endpoint.
"""

from chemsmart.agent.local.adapter import plan_to_synthesis_result
from chemsmart.agent.local.postprocessor import postprocess

__all__ = ["plan_to_synthesis_result", "postprocess"]
