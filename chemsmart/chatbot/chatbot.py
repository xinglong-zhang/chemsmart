from .executor import CommandExecutor
from .llm import LLM
from .llm_trainer import ChemSmartTrainer
from .logger import log


class Chatbot:
    def __init__(self):
        self.llm = LLM()
        self.executor = CommandExecutor()

    def chats(self):
        log("Chatbot started. Type 'exit' to quit.")
        while True:
            user_input = input("You: ")
            if user_input.lower() == "exit":
                log("Chatbot session ended.")
                break

            response = self.llm.generate_response(
                f"User command: {user_input}\nOutput:"
            )
            print(f"LLM: {response}")

            if response.startswith("EXECUTE:"):
                command = response.replace("EXECUTE:", "").strip()
                output = self.executor.execute(command)
                print(f"Command Output:\n{output}")

    def train(self):
        """Trains the LLM model."""
        trainer = ChemSmartTrainer()
        trainer.preprocess_data()
        trainer.train()


if __name__ == "__main__":
    Chatbot().chats()
