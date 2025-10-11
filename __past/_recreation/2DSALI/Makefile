# (1)コンパイラ
CC  = g++
# (2)コンパイルオプション
CFLAGS    = -Wall -Wextra -Wshadow -ftrapv -fstack-protector-all -Os
# (3)実行ファイル名
TARGET  = 2D_SALI
# (4)ディレクトリの定義
SRCDIR  = ../../src
INCDIR  = ../../inc
OBJDIR  = ../../obj
BINDIR  = ../../bin

# (5)コンパイル対象のソースコード
SRCS    = $(SRCDIR)/2D_SALI.cpp 
SRCS   += $(INCDIR)/calc.c $(INCDIR)/pcrtbp.cpp

# (6)オブジェクトファイル名
OBJS    = $(OBJDIR)/2D_SALI.o 
OBJS   += $(OBJDIR)/calc.o $(OBJDIR)/pcrtbp.o

# (7)インクルードディレクトリ
INCFLAGS = -I$(INCDIR)

# (8)VPATHの設定
VPATH = $(SRCDIR) $(INCDIR) $(OBJDIR)

# (9)ターゲットファイル生成
$(BINDIR)/$(TARGET): $(OBJS)
	$(CC) -o $@ $^

# (10)オブジェクトファイル生成
$(OBJDIR)/2D_SALI.o: 2D_SALI.cpp
	$(CC) $(CFLAGS) $(INCFLAGS) -c $< -o $@

$(OBJDIR)/calc.o: calc.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c $< -o $@

$(OBJDIR)/pcrtbp.o: pcrtbp.cpp
	$(CC) $(CFLAGS) $(INCFLAGS) -c $< -o $@

# (11)"make all"で make cleanとmakeを同時に実施。
all: clean $(BINDIR)/$(TARGET)

# (12).oファイル、実行ファイル、.dファイルを削除
clean:
	del /Q $(OBJS) $(BINDIR)/$(TARGET) *.d

# (13) .PHONYターゲットの追加
.PHONY: all clean
