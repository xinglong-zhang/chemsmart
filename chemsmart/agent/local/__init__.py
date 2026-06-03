"""Local LoRA model adapter for chemsmart agent.

Wraps the Smilesjs/chemsmart-qwen2.5-7b-lora model so that natural-language
queries can be served through the same ``chemsmart agent ask`` flow used for
OpenAI/Anthropic providers, without exposing remote API keys.

Submodules:
    loader        Load Qwen2.5-7B-Instruct base + LoRA adapter with 4-bit NF4.
    generator     Greedy decode wrapper that returns the parsed planner JSON.
    postprocessor Token-level repairs for V4 planner output.
    adapter       Convert atomic planner JSON to SynthesisSession schema.
    server        OpenAI-compatible /v1/chat/completions FastAPI endpoint.
"""

from chemsmart.agent.local.adapter import plan_to_synthesis_result
from chemsmart.agent.local.postprocessor import postprocess

__all__ = ["plan_to_synthesis_result", "postprocess"]
