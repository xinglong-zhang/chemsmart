"""Tests for the ChemSmart chatbot module."""


from chemsmart.chatbot.chatbot import Chatbot
from chemsmart.chatbot.llm import LLM


class TestLLM:
    """Tests for the LLM class."""

    def test_llm_init_offline_mode(self):
        """Test LLM initializes in offline mode without API key."""
        llm = LLM()
        assert llm.offline_mode is True

    def test_llm_loads_training_data(self):
        """Test LLM loads training data from commands_dataset.json."""
        llm = LLM()
        assert len(llm.training_data) > 0

    def test_llm_generate_response_offline(self):
        """Test LLM generates response in offline mode."""
        llm = LLM()
        response = llm.generate_response("Run optimization on test.xyz")
        assert "chemsmart" in response
        assert "opt" in response

    def test_llm_extract_file_from_input(self):
        """Test LLM extracts file name from user input."""
        llm = LLM()
        response = llm.generate_response(
            "Submit job for molecule.xyz on server shared"
        )
        assert "molecule.xyz" in response

    def test_llm_detect_job_type_ts(self):
        """Test LLM detects transition state job type."""
        llm = LLM()
        response = llm.generate_response(
            "Run transition state search on ts.com"
        )
        assert " ts" in response

    def test_llm_detect_job_type_sp(self):
        """Test LLM detects single point job type."""
        llm = LLM()
        response = llm.generate_response(
            "Submit single point calculation on test.log"
        )
        assert " sp" in response

    def test_llm_detect_job_type_irc(self):
        """Test LLM detects IRC job type."""
        llm = LLM()
        response = llm.generate_response("Run irc job on ts.log")
        assert " irc" in response

    def test_llm_detect_job_type_scan(self):
        """Test LLM detects scan job type."""
        llm = LLM()
        response = llm.generate_response("Run pes scan job on molecule.xyz")
        assert " scan" in response

    def test_llm_detect_server(self):
        """Test LLM detects server name from input."""
        llm = LLM()
        response = llm.generate_response(
            "Run opt job on test.com with server myserver"
        )
        assert "-s myserver" in response

    def test_llm_detect_project(self):
        """Test LLM detects project name from input."""
        llm = LLM()
        response = llm.generate_response(
            "Run opt job on test.com with project myproject"
        )
        assert "-p myproject" in response

    def test_llm_chat_method(self):
        """Test LLM chat method with message history."""
        llm = LLM()
        messages = [
            {"role": "user", "content": "Run optimization on test.xyz"}
        ]
        response = llm.chat(messages)
        assert "chemsmart" in response


class TestChatbot:
    """Tests for the Chatbot class."""

    def test_chatbot_init(self):
        """Test Chatbot initializes correctly."""
        chatbot = Chatbot()
        assert chatbot.llm is not None
        assert chatbot.executor is not None
        assert chatbot.auto_execute is False

    def test_chatbot_init_auto_execute(self):
        """Test Chatbot initializes with auto_execute option."""
        chatbot = Chatbot(auto_execute=True)
        assert chatbot.auto_execute is True

    def test_chatbot_is_chemsmart_command(self):
        """Test Chatbot identifies chemsmart commands."""
        chatbot = Chatbot()
        assert chatbot._is_chemsmart_command("chemsmart sub") is True
        assert chatbot._is_chemsmart_command("CHEMSMART sub") is True
        assert chatbot._is_chemsmart_command("other command") is False

    def test_chatbot_extract_command(self):
        """Test Chatbot extracts command from response."""
        chatbot = Chatbot()

        # Test direct command
        cmd = chatbot._extract_command("chemsmart sub -s shared gaussian opt")
        assert cmd == "chemsmart sub -s shared gaussian opt"

        # Test EXECUTE: prefix
        cmd = chatbot._extract_command(
            "EXECUTE: chemsmart sub -s shared gaussian opt"
        )
        assert cmd == "chemsmart sub -s shared gaussian opt"

        # Test command in backticks
        cmd = chatbot._extract_command(
            "Use `chemsmart sub -s shared gaussian opt`"
        )
        assert cmd == "chemsmart sub -s shared gaussian opt"

    def test_chatbot_generate_command(self):
        """Test Chatbot generates command from description."""
        chatbot = Chatbot()
        cmd = chatbot.generate_command("Run optimization on test.xyz")
        assert cmd is not None
        assert "chemsmart" in cmd

    def test_chatbot_generate_command_with_details(self):
        """Test Chatbot generates command with specific details."""
        chatbot = Chatbot()
        cmd = chatbot.generate_command(
            "Submit ts search on ts.com with server shared and project test"
        )
        assert cmd is not None
        assert "chemsmart" in cmd
        assert "-s shared" in cmd
        assert "-p test" in cmd
        assert " ts" in cmd

    def test_chatbot_process_message(self):
        """Test Chatbot processes a message and returns response."""
        chatbot = Chatbot()
        response, command, output = chatbot.process_message(
            "Run optimization on test.xyz"
        )
        assert response is not None
        # Command may or may not be extracted depending on response format
        # but conversation history should be updated
        assert len(chatbot.conversation_history) == 2

    def test_chatbot_conversation_history(self):
        """Test Chatbot maintains conversation history."""
        chatbot = Chatbot()
        chatbot.process_message("Hello")
        chatbot.process_message("Run optimization")
        assert len(chatbot.conversation_history) == 4
        assert chatbot.conversation_history[0]["role"] == "user"
        assert chatbot.conversation_history[1]["role"] == "assistant"


class TestCommandExecutor:
    """Tests for the CommandExecutor class."""

    def test_executor_execute_command(self):
        """Test CommandExecutor executes a command."""
        from chemsmart.chatbot.executor import CommandExecutor

        executor = CommandExecutor()
        result = executor.execute("echo hello")
        assert "hello" in result

    def test_executor_execute_invalid_command(self):
        """Test CommandExecutor handles invalid command."""
        from chemsmart.chatbot.executor import CommandExecutor

        executor = CommandExecutor()
        result = executor.execute("nonexistent_command_xyz123")
        assert "Error" in result or "not found" in result.lower()
