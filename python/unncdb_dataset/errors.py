class StructureRejected(ValueError):
    """A known data-quality condition, safe to record and continue past."""
    def __init__(self, code, message):
        super().__init__(message)
        self.code = code


class InsufficientDataError(RuntimeError):
    """The requested size/balance cannot be achieved under the selected conditions."""
    def __init__(self, message, report=None):
        super().__init__(message)
        self.report = report or {}
        self.partial_dataset = None
