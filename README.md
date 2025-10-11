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