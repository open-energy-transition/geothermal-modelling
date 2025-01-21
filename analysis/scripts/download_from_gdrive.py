# coding=utf-8# -*- coding: utf-8 -*-
# # SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
# #
# # SPDX-License-Identifier: AGPL-3.0-or-later
#
# # -*- coding: utf-8 -*-

import gdown
import pathlib
import requests
from datetime import datetime
from dateutil.relativedelta import relativedelta

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_usa import mock_snakemake

        snakemake = mock_snakemake("download_from_gdrive")

    url = snakemake.params["gdrive_url"]
    default_path = pathlib.Path(__file__).parent.parent.parent
    cookies_path = pathlib.Path(pathlib.Path.home(), snakemake.params["cookies_path"])
    pathlib.Path(cookies_path).mkdir(parents=True, exist_ok=True)
    cookies_file_path = pathlib.Path(cookies_path, "cookies.txt")
    cookies_file_path.touch(exist_ok=True)
    cookies_file = open(cookies_file_path, "w")
    download_path = pathlib.Path(default_path, snakemake.params["output_directory"])

    # get the cookies.txt file

    # --> send HTTP requests
    response = requests.get(url)

    # --> get response about Cookies
    cookies = response.cookies

    # --> store the cookies to cookies.txt
    cookies_file.write("# Netscape HTTP Cookie File \n")
    cookies_file.write("# http://curl.haxx.se/rfc/cookie_spec.html \n")
    cookies_file.write("# This is a generated file!  Do not edit. \n")
    cookies_file.write(" \n")

    domain = ".google.com"
    domain_specified = "TRUE"
    path = "/"
    secure = "TRUE"

    # add five months to today. Convert the time to epoch time. Cast it to integer and string
    expires = str(
        int(
            (
                datetime.today() + relativedelta(months=snakemake.params.delta_months)
            ).timestamp()
        )
    )

    cookie_list = [domain, domain_specified, path, secure, expires]

    for cookie in cookies:
        cookie_list.append(cookie.name)
        cookie_list.append(cookie.value)

    cookies_file.write("\t".join(cookie_list))

    cookies_file.close()

    # Download the data
    try:
        gdown.download_folder(
            url,
            output=str(download_path),
            quiet=False,
            resume=True,
            use_cookies=True,
        )
    except Exception as e:
        print("Error", e)
        pass
