# coding=utf-8# -*- coding: utf-8 -*-
# # SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
# #
# # SPDX-License-Identifier: AGPL-3.0-or-later
#
# # -*- coding: utf-8 -*-

import gdown
import pathlib
import requests

url = snakemake.params["gdrive_url"]
default_path = pathlib.Path(__file__).parent.parent.parent
cookies_path = pathlib.Path(default_path, snakemake.params["cookies_path"])
pathlib.Path(cookies_path).mkdir(parents=True, exist_ok=True)
cookies_file_path = pathlib.Path(
    cookies_path, "cookies.txt"
)
cookies_file_path.touch(exist_ok=True)
download_path = pathlib.Path(default_path, snakemake.params["output_path"])

# get the cookies.txt file

# --> send HTTP requests
response = requests.get(url)

# --> get response about Cookies
cookies = response.cookies

# --> store the cookies to cookies.txt
with open(cookies_file_path, 'w') as file:
    for cookie in cookies:
        file.write(f"{cookie.name}={cookie.value}\n")

# Download the data
try:
    gdown.download_folder(
        url, output=str(download_path), quiet=False, resume=True, use_cookies=True
    )
except Exception as e:
    print("Error", e)
    pass
