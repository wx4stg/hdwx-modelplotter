#!/usr/bin/env python3
# Purges no-longer-needed files from modelplotter
# Created on 19 December 2021 by Sam Gardner <stgardner4@tamu.edu>

from datetime import datetime as dt, timedelta
from os import path, walk, remove

def writeToStatus(stringToWrite):
    print(stringToWrite)
    stringToWrite = stringToWrite+"\n"
    with open(path.join(basePath, "status.txt"), "a") as statw:
        statw.write(stringToWrite)
        statw.close()

if __name__ == "__main__":
    now = dt.now()
    basePath = path.dirname(path.abspath(__file__))
    now = dt.now()
    if path.exists(path.join(basePath, "lastCleanDT.txt")):
        readLastCleanFile = open(path.join(basePath, "lastCleanDT.txt"), "r")
        lastCleanTime = dt.strptime(readLastCleanFile.read(), "%Y%m%d%H%M")
        readLastCleanFile.close()
    else:
        lastCleanTime = dt.utcnow()
        writeLastCleanFile = open(path.join(basePath, "lastCleanDT.txt"), "w")
        writeLastCleanFile.write(lastCleanTime.strftime("%Y%m%d%H%M"))
        writeLastCleanFile.close()
    if lastCleanTime < dt.utcnow() - timedelta(hours=2):
        writeToStatus("Cleaning up...")
        modelDataPath = path.join(basePath, "modelData")
        outputPath = path.join(basePath, "output")
        if path.exists(modelDataPath):
            for root, dirs, files in walk(modelDataPath):
                for name in files:
                    filepath = path.join(path.join(basePath, root), name)
                    createTime = dt.fromtimestamp(path.getmtime(filepath))
                    if createTime < now - timedelta(hours=3):
                        remove(filepath)
                        writeToStatus(filepath+" deleted.")
                    if filepath.endswith(".idx"):
                        remove(filepath)
                        writeToStatus(filepath+" deleted.")
        if path.exists(outputPath):
            for root, dirs, files in walk(outputPath):
                for name in files:
                    filepath = path.join(path.join(basePath, root), name)
                    if filepath.endswith(".json"):
                        deleteAfter = timedelta(days=3)
                    else:
                        deleteAfter = timedelta(minutes=20)
                    createTime = dt.fromtimestamp(path.getmtime(filepath))
                    if createTime < now - deleteAfter:
                        remove(filepath)
                        writeToStatus(filepath+" deleted.")
        remove(path.join(basePath, "status.txt"))
        remove(path.join(basePath, "lastCleanDT.txt"))