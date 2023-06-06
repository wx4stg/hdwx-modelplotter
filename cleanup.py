#!/usr/bin/env python3
# Purges no-longer-needed files from modelplotter
# Created on 19 December 2021 by Sam Gardner <stgardner4@tamu.edu>

from datetime import datetime as dt, timedelta
from os import path, walk, remove

if __name__ == "__main__":
    basePath = path.dirname(path.abspath(__file__))
    now = dt.now()
    modelDataPath = path.join(basePath, "modelData")
    outputPath = path.join(basePath, "output")
    if path.exists(modelDataPath):
        for root, dirs, files in walk(modelDataPath):
            for name in files:
                filepath = path.join(basePath, root, name)
                createTime = dt.fromtimestamp(path.getmtime(filepath))
                if createTime < now - timedelta(hours=3):
                    remove(filepath)
                if filepath.endswith(".idx"):
                    remove(filepath)
    if path.exists(outputPath):
        for root, dirs, files in walk(outputPath):
            for name in files:
                filepath = path.join(basePath, root, name)
                if filepath.endswith(".json"):
                    deleteAfter = timedelta(days=3)
                else:
                    deleteAfter = timedelta(minutes=20)
                createTime = dt.fromtimestamp(path.getmtime(filepath))
                if createTime < now - deleteAfter:
                    remove(filepath)