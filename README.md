# Meta-material Damping

A Python application for meta-material damping analysis and simulation.

## Project Structure

```
Meta-material-damping/
├── Application/
│   ├── src/                 # Source code
│   │   └── __init__.py
│   └── test/                # Test files
│       ├── __init__.py
│       └── test_example.py
├── Ref_Doc/                 # Reference documents
├── venv/                    # Virtual environment
├── pyproject.toml           # Project configuration
├── requirements.txt         # Dependencies
├── README.md                # This file
└── .gitignore               # Git ignore rules
```

## Setup

### 1. Activate the Virtual Environment

**Windows (PowerShell):**
```powershell
.\venv\Scripts\Activate.ps1
```

**Windows (Command Prompt):**
```cmd
venv\Scripts\activate.bat
```

### 2. Install Dependencies

```bash
pip install -r requirements.txt
```

Or install in development mode:
```bash
pip install -e ".[dev,jupyter]"
```

## Usage

```python
from src import __version__

print(f"Version: {__version__}")
```

## Running Tests

```bash
pytest
```

With coverage:
```bash
pytest --cov=Application/src
```

## Code Quality

Format code:
```bash
black Application/
```

Lint code:
```bash
flake8 Application/
```

Type checking:
```bash
mypy Application/src/
```

## License

MIT License
