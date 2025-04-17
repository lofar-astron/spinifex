from __future__ import annotations

from importlib import resources
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from aiohttp import web


@pytest.fixture
def dummy_bytes() -> bytes:
    return b"Example content for testing."


@pytest.fixture
def test_data_path() -> Path:
    with resources.as_file(resources.files("spinifex.data.tests")) as datapath:
        return Path(datapath)


@pytest.fixture
async def mock_http_server(tmp_path, dummy_bytes, test_data_path: Path):
    tmp_path = Path(tmp_path)
    file_path = tmp_path / "test_file.dat"
    file_path.write_bytes(dummy_bytes)

    file_paths = {file_path.name: file_path}

    test_data_file = test_data_path.glob("*")
    for file in test_data_file:
        if file.is_file():
            file_paths[file.name] = file

    async def generic_handler(request):
        filename = request.match_info["filename"]
        if filename in file_paths:
            return web.FileResponse(path=file_paths[filename])
        return web.Response(status=404, text="File not found")

    app = web.Application()
    app.router.add_get("/{filename}", generic_handler)

    runner = web.AppRunner(app)
    await runner.setup()
    site = web.TCPSite(runner, "localhost", 0)
    await site.start()

    port = site._server.sockets[0].getsockname()[1]
    url = f"http://localhost:{port}"

    yield url

    await runner.cleanup()


@pytest.fixture
def mock_ftp_server(tmp_path, dummy_bytes):
    with patch("spinifex.download.FTP") as MockFTP:
        tmp_path = Path(tmp_path)
        file_path = tmp_path / "test_file.dat"
        file_path.write_bytes(dummy_bytes)
        file_paths = {file_path.name: file_path}
        with resources.as_file(
            resources.files("spinifex.data.tests")
        ) as test_data_path:
            test_data_file = test_data_path.glob("*")
            for file in test_data_file:
                if file.is_file():
                    file_paths[file.name] = file

        ftp_instance = MagicMock()
        MockFTP.return_value = ftp_instance

        ftp_instance.login.return_value = None
        ftp_instance.cwd.return_value = None
        ftp_instance.quit.return_value = None

        def retrbinary_side_effect(cmd, callback):
            # Extract filename from "RETR filename"
            _, filename = cmd.split()
            file_path = file_paths.get(filename)
            if file_path is not None:
                callback(file_path.read_bytes())
            else:
                msg = f"No such file: {filename}"
                raise FileNotFoundError(msg)

        ftp_instance.retrbinary.side_effect = retrbinary_side_effect

        yield ftp_instance
