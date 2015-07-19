#!/bin/bash

sccBmc -c status | tee >(grep VCC0 > voltage_level_0.txt) | tee >(grep VCC1 > voltage_level_1.txt)| tee >(grep VCC3>voltage_level_3.txt)|tee >(grep VCC4>voltage_level_4.txt) | tee >(grep VCC5>voltage_level_5.txt)| tee >(grep VCC7>voltage_level_7.txt)| tee>last_voltage
