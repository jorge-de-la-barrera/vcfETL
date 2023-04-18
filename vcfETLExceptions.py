class vcfETLExceptions(Exception):
    """Base class for exceptions in the vcfETL class."""

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class NonNormalizedVCFError(vcfETLExceptions):
    """Exception raised if the VCF file is not normalized.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        super().__init__(message)

