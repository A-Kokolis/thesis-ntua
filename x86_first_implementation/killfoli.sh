#!/bin/bash 

kill $(ps aux | grep '[c]kpt_infoli/infoli' | awk '{print $2}')
