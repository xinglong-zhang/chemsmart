import json

import torch
from datasets import Dataset
from transformers import (
    AutoModelForCausalLM,
    AutoTokenizer,
    DataCollatorForSeq2Seq,
    Trainer,
    TrainingArguments,
    pipeline,
)


class ChemSmartTrainer:
    def __init__(
        self,
        model_name="deepseek-ai/DeepSeek-R1-Distill-Qwen-32B",
        dataset_path="chemsmart_data.json",
    ):
        """
        Initializes the ChemSmartTrainer class using a remote Hugging Face model.

        :param model_name: Hugging Face model name (e.g., "deepseek-ai/DeepSeek-R1-Distill-Qwen-32B").
        :param dataset_path: Path to the training dataset in JSON format.
        """
        self.model_name = model_name
        self.dataset_path = dataset_path
        self.pipeline = None
        self.model = None
        self.tokenizer = None
        self.tokenized_datasets = None

        # Load dataset
        self._load_dataset()

        # Load tokenizer & model
        self._load_pipeline()  # Loads a pipeline for inference
        self._load_model()  # Loads a model for training

    def _load_dataset(self):
        """Loads the dataset from JSON and converts it to Hugging Face Dataset format."""
        with open(self.dataset_path, "r") as f:
            data = json.load(f)

        self.dataset = Dataset.from_dict(
            {
                "input": [d["input"] for d in data],
                "output": [d["output"] for d in data],
            }
        )

    def _load_pipeline(self):
        """Loads a Hugging Face pipeline for inference (zero-shot command generation)."""
        self.pipeline = pipeline("text-generation", model=self.model_name)

    def _load_model(self):
        """Loads the tokenizer and model for training."""
        self.tokenizer = AutoTokenizer.from_pretrained(self.model_name)
        self.model = AutoModelForCausalLM.from_pretrained(
            self.model_name, torch_dtype=torch.float16, device_map="auto"
        )

    def _tokenize_function(self, examples):
        """Tokenizes input-output pairs for training."""
        return self.tokenizer(
            examples["input"],
            text_target=examples["output"],
            padding="max_length",
            truncation=True,
        )

    def preprocess_data(self):
        """Tokenizes the dataset and stores it in self.tokenized_datasets."""
        self.tokenized_datasets = self.dataset.map(
            self._tokenize_function, batched=True
        )

    def train(self, output_dir="./chemsmart-llm", num_epochs=5, batch_size=4):
        """
        Fine-tunes the model using the prepared dataset.

        :param output_dir: Directory to save the trained model.
        :param num_epochs: Number of training epochs.
        :param batch_size: Batch size for training.
        """
        if self.tokenized_datasets is None:
            raise ValueError(
                "Dataset not tokenized. Call preprocess_data() before training."
            )

        training_args = TrainingArguments(
            output_dir=output_dir,
            per_device_train_batch_size=batch_size,
            per_device_eval_batch_size=batch_size,
            num_train_epochs=num_epochs,
            logging_dir="./logs",
            save_total_limit=2,
            save_strategy="epoch",
            evaluation_strategy="epoch",
            fp16=True,
        )

        trainer = Trainer(
            model=self.model,
            args=training_args,
            train_dataset=self.tokenized_datasets,
            eval_dataset=self.tokenized_datasets,
            tokenizer=self.tokenizer,
            data_collator=DataCollatorForSeq2Seq(self.tokenizer),
        )

        trainer.train()
        self.save_model(output_dir)

    def save_model(self, save_path="./chemsmart_model"):
        """Saves the fine-tuned model and tokenizer."""
        self.model.save_pretrained(save_path)
        self.tokenizer.save_pretrained(save_path)

    def generate_command(self, user_input):
        """
        Uses the pipeline to generate a chemistry job submission command from user input.

        :param user_input: User's text instruction.
        :return: Model-generated job submission command.
        """
        response = self.pipeline(
            user_input, max_length=100, num_return_sequences=1
        )
        return response[0]["generated_text"]


if __name__ == "__main__":
    trainer = ChemSmartTrainer()

    # Example inference (interactive command generation)
    user_input = "Please submit a DFT optimization job using Gaussian with project settings toluene for file input.xyz."
    generated_command = trainer.generate_command(user_input)
    print("\nGenerated Command:", generated_command)

    # Uncomment to train
    # trainer.preprocess_data()
    # trainer.train()
