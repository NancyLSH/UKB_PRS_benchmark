#!/bin/bash

SRC="/mnt/sdc1/UKB/Genotype/Chinese/"
DST="/data1/jiapl_group/lishuhua/UKB/Genotype/Chinese/"
LOGDIR="/data1/jiapl_group/lishuhua/UKB/logs/"

mkdir -p "$LOGDIR"

LOG="$LOGDIR/rsync_copy_3_$(date +%F).log"
MD5SRC="$LOGDIR/md5_source_3_$(date +%F).txt"
MD5DST="$LOGDIR/md5_target_3_$(date +%F).txt"


echo "Start time: $(date)" >> $LOG
echo "Source: $SRC" >> $LOG
echo "Destination: $DST" >> $LOG
echo "Log file: $LOG" >> $LOG
echo "MD5 Source file: $MD5SRC" >> $LOG
echo "MD5 Target file: $MD5DST" >> $LOG
echo ""

echo "Copying files from $SRC to $DST..." | tee -a $LOG
rsync -av --partial --progress "$SRC/" "$DST/" >> "$LOG" 2>&1

if [ $? -ne 0 ]; then
    echo "[!] rsync failed. Check the log file for details." | tee -a $LOG
    exit 1
fi

echo "Transfer completed successfully." | tee -a $LOG
echo "Calculating MD5 checksums..." | tee -a $LOG

SUCCESS=0
FAILURE=0
MISSING=0

> "$MD5SRC"
> "$MD5DST"

find "$SRC" -type f | while read -r src_file; do
    rel_path="${src_file#$SRC}"
    dst_file="$DST/$rel_path"

    if [ -f "$dst_file" ]; then
        src_md5=$(md5sum "$src_file" | awk '{print $1}')
        dst_md5=$(md5sum "$dst_file" | awk '{print $1}')
        # save md5 to a file
        echo "$src_md5 $rel_path" >> "$MD5SRC"
        echo "$dst_md5 $rel_path" >> "$MD5DST"
        if [ "$src_md5" == "$dst_md5" ]; then
            echo "MD5 match: $rel_path" | tee -a $LOG
            SUCCESS=$((SUCCESS + 1))
        else
            echo "[!] MD5 mismatch: $rel_path" | tee -a $LOG
            FAILURE=$((FAILURE + 1))
        fi
    else
        echo "[!] Missing file in destination: $rel_path" | tee -a $LOG
        MISSING=$((MISSING + 1))
    fi
done

echo "Summary:" | tee -a $LOG
echo "  Successful transfers: $SUCCESS" | tee -a $LOG
echo "  Failed transfers: $FAILURE" | tee -a $LOG
echo "  Missing files: $MISSING" | tee -a $LOG
echo "End time: $(date)" | tee -a $LOG