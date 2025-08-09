#!/bin/bash

# 设置新的内存要求(8GB)
NEW_MEMORY=8192

# 获取用户当前被挂起(H)的condor作业ID
HOLD_JOBS=$(condor_q -hold -format "%d." ClusterId -format "%d\n" ProcId)

# 检查是否有被挂起的作业
if [ -z "$HOLD_JOBS" ]; then
    echo "No held condor job found"
    exit 0
fi

echo "Found held job:"
echo "$HOLD_JOBS"

# 遍历每个被挂起的作业
for job in $HOLD_JOBS; do
    # 分割ClusterId和ProcId
    IFS='.' read -ra ID <<< "$job"
    cluster=${ID[0]}
    proc=${ID[1]}
    
    echo "processing: $cluster.$proc"
    
    condor_qedit $cluster.$proc RequestMemory "$NEW_MEMORY"

    condor_release $cluster.$proc
    
done

echo "All done"
