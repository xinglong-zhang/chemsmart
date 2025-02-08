import os

MODEL_PATH = os.path.join(
    os.path.dirname(__file__),
    "models/ggml-model-q4_0.gguf",
)  # local model

# MODEL_NAME = "mistral-7b-instruct-v0.1.Q4_K_M.gguf"
# MODEL_NAME = "Qwen/Qwen2-0.5B-Instruct"  # requires transformers>=4.37.0
# MODEL_NAME = "deepseek-ai/DeepSeek-R1-Distill-Qwen-32B"  # requires transformers>=4.37.0
MODEL_NAME = "TinyLlama/TinyLlama-1.1B-Chat-v1.0"
# MODEL_NAME = "microsoft/phi-2"

COMMANDS_DATASET = os.path.join(
    os.path.dirname(__file__),
    "training_data/commands_dataset.json",
)

CHEMSMART_MODEL = os.path.join(
    os.path.dirname(__file__),
    "models/chemsmart-llm",
)

print(MODEL_PATH)
print(MODEL_NAME)
print(COMMANDS_DATASET)
