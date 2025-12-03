awk -F, '$2 != 99 && $2 != 0' gongyou.csv > cleaned.csv     # 去除99和0
