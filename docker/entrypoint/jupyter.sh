#!/bin/bash
sage -n jupyter --no-browser --ip=$(grep `hostname` /etc/hosts | cut -f1) --port=8888
