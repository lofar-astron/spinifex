from __future__ import annotations

from importlib import resources
from pathlib import Path

import pytest
from aiohttp import web
from spinifex.download import download_file_http

EXAMPLE_CONTENT = b"Example content for testing."


@pytest.fixture
async def local_file_server(tmp_path):
    tmp_path = Path(tmp_path)
    file_path = tmp_path / "test_file.dat"
    file_path.write_bytes(EXAMPLE_CONTENT)

    file_paths = {file_path.name: file_path}

    with resources.as_file(resources.files("spinifex.data.tests")) as test_data_path:
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


@pytest.mark.allow_hosts(["127.0.0.1", "::1"])
@pytest.mark.asyncio
async def test_download_file_http(local_file_server, tmp_path):
    tmp_path = Path(tmp_path)
    file_name = "test_file.dat"
    url = f"{local_file_server}/{file_name}"

    output_file = tmp_path / "downloaded.dat"
    await download_file_http(url, output_file)

    # Assert file was written with expected content
    assert output_file.exists()
    assert output_file.read_bytes() == EXAMPLE_CONTENT


@pytest.mark.allow_hosts(["127.0.0.1", "::1"])
@pytest.mark.asyncio
async def test_download_file_http_test_data(local_file_server, tmp_path):
    tmp_path = Path(tmp_path)

    file_name = "casg0010.99i.Z"
    url = f"{local_file_server}/{file_name}"

    output_file = tmp_path / file_name
    await download_file_http(url, output_file)

    # Assert file was written with expected content
    assert output_file.exists()
    with resources.as_file(resources.files("spinifex.data.tests")) as test_data_path:
        test_data_file = test_data_path / file_name
        assert output_file.read_bytes() == test_data_file.read_bytes()
