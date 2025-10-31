#!/usr/bin/env python

import argparse
import logging
import os
import re

import boto3
from botocore import UNSIGNED
from botocore.config import Config

FIRE_ENDPOINT: str = "https://hl.fire.sdo.ebi.ac.uk"
PUBLIC_FTP_PATH: str = "ftp://ftp.sra.ebi.ac.uk/vol1/"
PRIVATE_FTP_PATH: str = "ftp://ftp.dcc-private.ebi.ac.uk/vol1/"
PUBLIC_BUCKET: str = "era-public"
PRIVATE_BUCKET: str = "era-private"


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def transform_ftp_to_s3(ftp_path: str) -> tuple[str, str]:
    """
    Transforms an FTP or HTTP path to a FIRE S3 object key, it also returns if it's public or private.
    :param ftp_path: The FTP path of the file to be transformed.
    :type ftp_path: str
    :return: A tuple containing the S3 object key and the corresponding bucket name.
    :rtype: tuple[str, str]
    :raises ValueError: If the FTP path does not match the expected format.
    """
    # Normalize protocol to ftp://
    normalized_path = re.sub(r"^(https?|ftp)://", "ftp://", ftp_path)
    if not normalized_path.startswith("ftp://"):
        normalized_path = f"ftp://{normalized_path}"

    # Check if public or private
    if normalized_path.startswith(PUBLIC_FTP_PATH):
        s3_key = normalized_path.replace(PUBLIC_FTP_PATH, "")
        logger.info(f"Detected a public file for FTP path: {normalized_path}")
        return s3_key, PUBLIC_BUCKET
    elif normalized_path.startswith(PRIVATE_FTP_PATH):
        s3_key = normalized_path.replace(PRIVATE_FTP_PATH, "")
        logger.info(f"Detected a private file for FTP path: {normalized_path}")
        return s3_key, PRIVATE_BUCKET
    else:
        raise ValueError(
            f"Invalid FTP path: {normalized_path}. Must start with {PUBLIC_FTP_PATH} or {PRIVATE_FTP_PATH}."
        )


def download_file_from_fire(
    s3_key: str,
    bucket: str,
    outdir: str,
    access_key: str | None = None,
    secret_key: str | None = None,
) -> None:
    """
    Downloads an individual file from FIRE S3 using its object key.

    :param s3_key: The S3 object key of the file to download.
    :type s3_key: str
    :param bucket: The name of the S3 bucket.
    :type bucket: str
    :param outdir: The local directory to save the downloaded file.
    :type outdir: str
    :param access_key: The access key for private S3 buckets (optional for public files).
    :type access_key: str | None
    :param secret_key: The secret key for private S3 buckets (optional for public files).
    :type secret_key: str | None
    :return: None
    :rtype: None
    :raises ValueError: If credentials are missing for private files.
    :raises Exception: For other download errors.
    """
    s3_args = {"endpoint_url": FIRE_ENDPOINT}
    if bucket == PRIVATE_BUCKET:
        if not access_key or not secret_key:
            logger.error("Missing credentials for private files.")
            raise ValueError(
                "Access key and secret key are required for private files."
            )
        s3_args.update(
            aws_access_key_id=access_key,
            aws_secret_access_key=secret_key,
        )
    else:
        # Public bucket configuration with unsigned requests
        s3_args.update({"config": Config(signature_version=UNSIGNED)})

    s3 = boto3.client("s3", **s3_args)

    os.makedirs(outdir, exist_ok=True)
    local_file_path = os.path.join(outdir, os.path.basename(s3_key))

    try:
        logger.info(
            f"Downloading {s3_key} from S3 bucket {bucket} to {local_file_path}..."
        )
        s3.download_file(bucket, s3_key, local_file_path)
        logger.info(f"File successfully downloaded to: {local_file_path}")
    except Exception as e:
        logger.error(f"Error downloading file from S3: {e}")
        raise


def download_files(
    ftp_paths: list[str],
    outdir: str,
    access_key: str | None,
    secret_key: str | None,
) -> None:
    """
    Downloads multiple files from their FTP paths.

    :param ftp_paths: List of FTP paths to download.
    :type ftp_paths: list[str]
    :param outdir: Directory to save the downloaded files.
    :type outdir: str
    :param access_key: Access key for private files.
    :type access_key: str | None
    :param secret_key: Secret key for private files.
    :type secret_key: str | None
    """
    for ftp_path in ftp_paths:
        s3_key, bucket = transform_ftp_to_s3(ftp_path)
        download_file_from_fire(s3_key, bucket, outdir, access_key, secret_key)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download multiple files from FTP paths via FIRE S3 (supports public and private files)."
    )
    parser.add_argument(
        "--ftp-paths",
        nargs="+",
        required=True,
        help=f"Space-separated list of FTP paths to download (e.g., {PUBLIC_FTP_PATH}.../file1 {PRIVATE_FTP_PATH}.../file2).",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Local destination directory for the downloaded files.",
    )
    parser.add_argument(
        "--access-key",
        required=False,
        help="S3 access key (required for private files).",
    )
    parser.add_argument(
        "--secret-key",
        required=False,
        help="S3 secret key (required for private files).",
    )
    args = parser.parse_args()

    logger.info("Starting the file download process...")
    download_files(args.ftp_paths, args.outdir, args.access_key, args.secret_key)
    logger.info("All files have been processed.")


if __name__ == "__main__":
    main()
