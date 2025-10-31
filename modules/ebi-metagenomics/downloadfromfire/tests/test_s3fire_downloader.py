#!/usr/bin/env python

import sys
from pathlib import Path

import pytest

# Add the script directory to the path
script_dir = Path(__file__).parent.parent / "resources" / "usr" / "bin"
sys.path.insert(0, str(script_dir))

# Import the module to test
from s3fire_downloader import (  # noqa: E402
    PRIVATE_BUCKET,
    PUBLIC_BUCKET,
    transform_ftp_to_s3,
)


class TestTransformFtpToS3:
    """Test cases for the transform_ftp_to_s3 function."""

    @pytest.mark.parametrize(
        "ftp_path,expected_key,expected_bucket",
        [
            # Public FTP path
            (
                "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR123456_1.fastq.gz",
                "fastq/ERR123/ERR123456_1.fastq.gz",
                PUBLIC_BUCKET,
            ),
            # Private FTP path
            (
                "ftp://ftp.dcc-private.ebi.ac.uk/vol1/fastq/ERR456/ERR456789_1.fastq.gz",
                "fastq/ERR456/ERR456789_1.fastq.gz",
                PRIVATE_BUCKET,
            ),
            # HTTP protocol normalization
            (
                "http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR123456_1.fastq.gz",
                "fastq/ERR123/ERR123456_1.fastq.gz",
                PUBLIC_BUCKET,
            ),
            # HTTPS protocol normalization
            (
                "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR123456_1.fastq.gz",
                "fastq/ERR123/ERR123456_1.fastq.gz",
                PUBLIC_BUCKET,
            ),
            # Path without protocol
            (
                "ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR123456_1.fastq.gz",
                "fastq/ERR123/ERR123456_1.fastq.gz",
                PUBLIC_BUCKET,
            ),
            # Nested directory structure (public)
            (
                "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR123/ERR123456/file.fastq.gz",
                "run/ERR123/ERR123456/file.fastq.gz",
                PUBLIC_BUCKET,
            ),
            # Nested directory structure (private)
            (
                "ftp://ftp.dcc-private.ebi.ac.uk/vol1/run/ERR456/ERR456789/file.fastq.gz",
                "run/ERR456/ERR456789/file.fastq.gz",
                PRIVATE_BUCKET,
            ),
            # Special characters in filename
            (
                "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/file_name-v1.2.fastq.gz",
                "fastq/ERR123/file_name-v1.2.fastq.gz",
                PUBLIC_BUCKET,
            ),
        ],
    )
    def test_valid_transformations(self, ftp_path, expected_key, expected_bucket):
        """Test valid FTP path transformations with various formats."""
        s3_key, bucket = transform_ftp_to_s3(ftp_path)

        assert s3_key == expected_key
        assert bucket == expected_bucket

    @pytest.mark.parametrize(
        "invalid_path,expected_error_msg",
        [
            # Invalid server
            (
                "ftp://invalid.server.com/vol1/file.fastq.gz",
                "Invalid FTP path",
            ),
            # Missing vol1
            (
                "ftp://ftp.sra.ebi.ac.uk/fastq/ERR123/ERR123456_1.fastq.gz",
                "Invalid FTP path",
            ),
        ],
    )
    def test_invalid_paths_raise_error(self, invalid_path, expected_error_msg):
        """Test that invalid FTP paths raise ValueError."""
        with pytest.raises(ValueError) as exc_info:
            transform_ftp_to_s3(invalid_path)

        assert expected_error_msg in str(exc_info.value)

    def test_return_type_is_tuple(self):
        """Test that the function returns a tuple of two strings."""
        ftp_path = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR123456_1.fastq.gz"

        result = transform_ftp_to_s3(ftp_path)

        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], str)
        assert isinstance(result[1], str)
