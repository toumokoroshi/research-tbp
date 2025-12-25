# tabata repo readme

## adgenda

1. ディレクトリ構成とファイルの役割
2. ビルド方法
3. シミュレータの使用方法

## main

### 1.ディレクトリ構成とファイルの役割

TBP/
├── .clang-format             # c/cppフォーマッタ設定　　コードの動作やビルドの挙動には関係ない
├── .clangd　                  # c/cpp言語サーバーclangdの設定　　コードの動作やビルドの挙動には関係ない
├── .gitignore　               # gitignore　　コードの動作やビルドの挙動には関係ない
├── .gitmessage.txt           # コミットメッセージテンプレート　　コードの動作やビルドの挙動には関係ない
├── compile_flags.txt　        # clangdの設定として必要なはず　　コードの動作やビルドの挙動には関係ない
├── CMakeLists.txt            # ルートのCMake設定
├── CMakePresets.json         # CMakeのビルドのconfig
├── README.md　                # README
├── setup.bat                 # インテルコンパイラ用の環境変数設定を行うバッチファイルを走らせた後に、宣言されているすべてのappをcmakeを使いビルドするためのバッチファイル
│
├── .vscode/
│   └── settings.json         # vscodeの設定　フォーマッタの設定ファイルパスを示してあるだけ　コードの動作やビルドの挙動には関係ない
├── .zed/
│   └── settings.json　        # zedというエディタの設定　コードの動作やビルドの挙動には関係ない
├── __past/
│
├── _bin/
│
├── apps/                     # シミュレーション本体(app)たちを格納する
│   ├── CMakeLists.txt        # appを宣言するCMake
│   ├── SALI2d/
│   │   ├── CMakeLists.txt    # appのビルドに必要なライブラリをリンクするなどappのビルドに必要な個別設定
│   │   └── main.cpp
│   ├── SALI3d/
│   │   ├── CMakeLists.txt    # appのビルドに必要なライブラリをリンクするなどappのビルドに必要な個別設定
│   │   └── main.cpp
│   ├── SALI3dV2/
│   │   ├── CMakeLists.txt　　　　# appのビルドに必要なライブラリをリンクするなどappのビルドに必要な個別設定
│   │   └── main.cpp
│   └── SSBI2SEBR/
│       ├── CMakeLists.txt　　　　# appのビルドに必要なライブラリをリンクするなどappのビルドに必要な個別設定
│       └── main.cpp
├── configs/                  # シミュレータが読む設定ファイル
│   ├── 3D_crtbp_SALI/
│   │    └──  .....txt
│   └── astro_param
│        └── astro_param.txt  # JPLのデータベースから引っ張ってきた太陽とか地球のデータ
├── data/                  # シミュレータが読む設定ファイル
│   ├── 3D_crtbp_SALI/
│   │
│   ├── Ephemeris/
│   │
│   ├── result/
│   │
│   └── SALI/
├── lib/                      # 共有ライブラリ
│   ├── CMakeLists.txt
│   ├── include/              # 公開ヘッダーファイル
│   │   ├── crtbp.h           # 円制限三体問題での計算用のクラス　実装が汚いのでもう使っていない　テンプレートを使っているためヘッダオンリー
│   │   ├── rtbp.h　　　　　　　　　　　　# 円制限三体問題での計算用のクラス　テンプレートを使っているためヘッダオンリー
│   │   ├── utils.h　　　　　　　　　　　# プログレスバーなど汎用の関数をまとめたファイル　テンプレートを使っているためヘッダオンリー
│   │   └── vector3d.h　　　　　　　　# 基本的な3次元ベクトル演算を行うためのベクトルクラス　ほんとはeigenとか使うべき　テンプレートを使っているためヘッダオンリー
│   └── src/
│       └── .gitkeep
├── scripts/
│   ├── requirements.txt
│   └── README.md
├── test/                     # いろんなテストやデバッグが入ってる
│   ├── .gitkeep
│   └── README.md
└── external/                 # 外部依存ライブラリ（オプション）
    └── .gitkeep

### 2.ビルド方法
https://dexall.co.jp/articles/?p=1968
### 3.シミュレータの使用方法
