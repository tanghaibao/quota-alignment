#!/bin/bash

cluster_utils.py --dag --precision 1 dag_file cluster_file
quota_align.py --merge --quota 1:1 cluster_file
