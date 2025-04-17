from __future__ import annotations

from pathlib import Path

import pytest
from spinifex.download import download_file_ftp, download_file_http


@pytest.mark.allow_hosts(["127.0.0.1", "::1"])
@pytest.mark.asyncio
async def test_download_file_ftp(mock_ftp_server, tmp_path, dummy_bytes):
    _ = mock_ftp_server
    tmp_path = Path(tmp_path)
    output_file = tmp_path / "ftp_test_file.dat"
    url = "ftp://mockserver/some/dir/test_file.dat"

    await download_file_ftp(url, output_file)

    assert output_file.exists()
    assert output_file.read_bytes() == dummy_bytes

    output_file.unlink(missing_ok=True)


@pytest.mark.allow_hosts(["127.0.0.1", "::1"])
@pytest.mark.asyncio
async def test_download_file_http(mock_http_server, tmp_path, dummy_bytes):
    tmp_path = Path(tmp_path)
    file_name = "test_file.dat"
    url = f"{mock_http_server}/{file_name}"

    output_file = tmp_path / "http_test_file.dat"
    await download_file_http(url, output_file)

    # Assert file was written with expected content
    assert output_file.exists()
    assert output_file.read_bytes() == dummy_bytes

    output_file.unlink(missing_ok=True)


@pytest.mark.allow_hosts(["127.0.0.1", "::1"])
@pytest.mark.asyncio
async def test_download_file_http_test_data(mock_http_server, tmp_path, test_data_path):
    tmp_path = Path(tmp_path)

    file_name = "casg0010.99i.Z"
    url = f"{mock_http_server}/{file_name}"

    output_file = tmp_path / file_name
    await download_file_http(url, output_file)

    # Assert file was written with expected content
    assert output_file.exists()
    test_data_file = test_data_path / file_name
    assert output_file.read_bytes() == test_data_file.read_bytes()

    output_file.unlink(missing_ok=True)
