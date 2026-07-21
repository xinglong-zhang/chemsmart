import subprocess
import time


def run_accumulation(model, provider, corpus_name, limit):
    print(f"[{time.strftime('%X')}] Running {model} on {corpus_name}...")
    cmd = [
        "python",
        "scripts/training/reasoning_accum.py",
        model,
        "--provider",
        provider,
        "--corpus",
        corpus_name,
        "--limit",
        str(limit),
    ]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running {corpus_name} with {model}: {e}")


def main():
    # 3 Rounds of execution
    for round_num in range(1, 4):
        print(f"=== Starting Round {round_num} ===")
        for batch_num in range(1, 9):
            corpus_name = f"batch_{batch_num}"
            # OpenAI
            run_accumulation("gpt-5.4-luna", "openai", corpus_name, 5)
            # DeepSeek
            run_accumulation("deepseek-v4-pro", "deepseek", corpus_name, 5)

        print(f"=== Round {round_num} Complete. Exporting === ")
        subprocess.run(
            ["python", "scripts/training/export_sft.py", "--include-runs"],
            check=True,
        )


if __name__ == "__main__":
    main()
