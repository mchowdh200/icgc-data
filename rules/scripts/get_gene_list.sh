#!/bin/env bash
# input is piped in
awk 'BEGIN{OFS="\t"}{print $24,$25,$26,$27,$11}' | sort | uniq
