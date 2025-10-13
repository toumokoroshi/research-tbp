# tabata repo readme

## adgenda

1. ディレクトリ構成とファイルの役割
2. ビルド方法
3. シミュレータの使用方法

## main

### 1.ディレクトリ構成とファイルの役割

project_root/
├── CMakeLists.txt                 # ルートのCMake設定
├── README.md
├── .gitignore
│
├── lib/                           # 共有ライブラリ
│   ├── CMakeLists.txt
│   ├── include/                   # 公開ヘッダーファイル
│   │   └── mylib/
│   │       ├── math_utils.h
│   │       ├── solver.h
│   │       └── io_utils.h
│   └── src/                       # ライブラリの実装
│       ├── math_utils.cpp
│       ├── solver.cpp
│       └── io_utils.cpp
│
├── apps/                          # 個別のmain関数を持つアプリケーション
│   ├── CMakeLists.txt
│   ├── simulator1/
│   │   ├── CMakeLists.txt
│   │   └── main.cpp
│   ├── simulator2/
│   │   ├── CMakeLists.txt
│   │   └── main.cpp
│   └── analyzer/
│       ├── CMakeLists.txt
│       └── main.cpp
│
├── configs/                       # 設定ファイル
│   ├── default.conf
│   ├── scenario1.conf
│   └── scenario2.conf
│
├── output/                        # シミュレーション出力（gitignore推奨）
│   ├── .gitkeep
│   └── README.md
│
├── tests/                         # テストコード（オプション）
│   ├── CMakeLists.txt
│   └── test_*.cpp
│
└── external/                      # 外部依存ライブラリ（オプション）
    └── .gitkeep

### 2.ビルド方法

### 3.シミュレータの使用方法

`simulator1`だけをビルドする方法を説明します。

## 方法1: CMakeのターゲット指定（推奨）

```bash
# ビルドディレクトリで
cd build

# simulator1だけをビルド
cmake --build . --target simulator1

# または makeコマンド（Makefileジェネレータの場合）
make simulator1
```

## 方法2: 初回ビルドから全て

```bash
# プロジェクトルートから
mkdir build
cd build

# CMake設定
cmake ..

# simulator1だけをビルド
cmake --build . --target simulator1

# 実行
./bin/simulator1 ../configs/scenario1.conf
```

## 方法3: 特定のアプリだけビルドするオプション

apps/CMakeLists.txtを修正して、ビルドするアプリを選択可能にする:

### apps/CMakeLists.txt（修正版）

```cmake
# ビルドオプションの追加
option(BUILD_SIMULATOR1 "Build simulator1" ON)
option(BUILD_SIMULATOR2 "Build simulator2" ON)
option(BUILD_ANALYZER "Build analyzer" ON)

# 条件付きでサブディレクトリを追加
if(BUILD_SIMULATOR1)
    add_subdirectory(simulator1)
endif()

if(BUILD_SIMULATOR2)
    add_subdirectory(simulator2)
endif()

if(BUILD_ANALYZER)
    add_subdirectory(analyzer)
endif()
```

### 使い方

```bash
# simulator1だけをビルド対象にする
cmake .. -DBUILD_SIMULATOR1=ON -DBUILD_SIMULATOR2=OFF -DBUILD_ANALYZER=OFF

# または
cmake .. -DBUILD_SIMULATOR2=OFF -DBUILD_ANALYZER=OFF

# ビルド
cmake --build .
```

## 方法4: 並列ビルドで高速化

```bash
# 4スレッドで並列ビルド
cmake --build . --target simulator1 -j4

# または
make simulator1 -j4
```

## すべてのターゲット確認

```bash
# ビルド可能なターゲット一覧を表示
cmake --build . --target help

# または
make help
```

出力例:
```
The following are some of the valid targets for this Makefile:
... all (the default if no target is provided)
... clean
... depend
... mylib
... simulator1
... simulator2
... analyzer
```

## よくあるワークフロー

### 開発中（simulator1のみ修正・ビルド）

```bash
cd build

# simulator1だけビルド（高速）
cmake --build . --target simulator1 -j4

# すぐ実行
./bin/simulator1 ../configs/scenario1.conf
```

### クリーンビルド

```bash
# simulator1だけクリーン＆ビルド
cd build
rm -rf CMakeFiles/simulator1.dir  # simulator1の中間ファイルを削除
cmake --build . --target simulator1
```

### 全体をクリーンビルド

```bash
# buildディレクトリごと削除して再ビルド
cd ..
rm -rf build
mkdir build && cd build
cmake ..
cmake --build . --target simulator1
```

## Intel oneAPIコンパイラを使う場合

```bash
# oneAPI環境設定
source /opt/intel/oneapi/setvars.sh

# ビルドディレクトリ作成
mkdir build && cd build

# Intelコンパイラ指定
cmake -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx ..

# simulator1だけビルド
cmake --build . --target simulator1 -j4

# 実行
./bin/simulator1 ../configs/scenario1.conf
```

## デバッグビルド vs リリースビルド

```bash
# デバッグビルド（最適化なし、デバッグ情報あり）
cmake -DCMAKE_BUILD_TYPE=Debug ..
cmake --build . --target simulator1

# リリースビルド（最適化あり）
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --target simulator1
```

## まとめ

**最もシンプルな方法:**
```bash
cd build
cmake --build . --target simulator1
./bin/simulator1 ../configs/scenario1.conf
```

これで`simulator1`だけが再ビルドされ、他のアプリ（simulator2, analyzer）には影響しません！

project_root/
├── CMakeLists.txt                 # ルートのCMake設定
├── README.md
├── .gitignore
│
├── lib/                           # 共有ライブラリ
│   ├── CMakeLists.txt
│   ├── include/                   # 公開ヘッダーファイル
│   │   └── mylib/
│   │       ├── lib1.h
│   │       ├── lib2.h # headeronly
│   │       └── lib3.h # headeronly
│   └── src/                       # ライブラリの実装
│       └── lib1.cpp
│
├── apps/                          # 個別のmain関数を持つアプリケーション
│   ├── CMakeLists.txt
│   ├── simulator1/
│   │   ├── CMakeLists.txt
│   │   └── main.cpp
│   └── simulator2/
│       ├── CMakeLists.txt
│       └── main.cpp
│
├── configs/                       # 設定ファイル
│   ├── default.conf
│   ├── scenario1.conf
│   └── scenario2.conf
│
├── output/                        # シミュレーション出力（gitignore推奨）
│   ├── .gitkeep
│   └── README.md
│
└── external/
    ├── eigen/          
    ├── boost/          

    # oneAPI環境のセットアップ
source /opt/intel/oneapi/setvars.sh

# ビルドディレクトリの作成と移動
mkdir -p build
cd build
    CC=icx CXX=icpx cmake ..
cmake --build . --target simulator1

