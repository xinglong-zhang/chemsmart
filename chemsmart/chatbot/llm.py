from llama_cpp import Llama

from . import MODEL_PATH


class LLM:
    def __init__(self, model_path=MODEL_PATH):
        self.model = Llama(
            model_path=model_path, n_ctx=2048, n_threads=4, verbose=False
        )

    def generate_response(self, prompt):
        response = self.model(prompt, max_tokens=256, stop=["\n"], echo=False)
        return response["choices"][0]["text"].strip()
