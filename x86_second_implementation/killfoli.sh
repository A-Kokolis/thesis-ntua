#!/bin/bash 

kill $(ps aux | grep '[d]vfs_thread/infoli' | awk '{print $2}')
