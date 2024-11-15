# coding=utf-8# -*- coding: utf-8 -*-
# # SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
# #
# # SPDX-License-Identifier: AGPL-3.0-or-later
#
# # -*- coding: utf-8 -*-

import gdown
import pathlib


url = "https://drive.google.com/drive/folders/1LwSoZDtnyUx5ki9SmBvdlGW3QwWjs4rA?usp=drive_link"
default_path = pathlib.Path(__file__).parent.parent.parent
download_path = str(pathlib.Path(default_path, "analysis", "gdrive_data"))

try:
    gdown.download_folder(url, output=download_path, quiet=False, resume=True, use_cookies=True)
except Exception as e:
    print("Error", e)
    pass
