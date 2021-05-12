#!/bin/env bash
grep RG $1 | head -1 | tr '\t' '\n' | grep SM | cut -d':' -f2

