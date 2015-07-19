#!/bin/bash 

kill $(ps aux | grep 'monitor3V3SCC' | awk '{print $2}')
