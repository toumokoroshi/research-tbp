# TBP - Three Body Problem Simulation

円制限三体問題（CR3BP）のシミュレーションおよび解析ツール群

## 目次

1. [クイックスタート](#クイックスタート)
2. [必要な環境](#必要な環境)
3. [手動セットアップ](#手動セットアップ)
4. [ディレクトリ構成](#ディレクトリ構成)
5. [シミュレータの使用方法](#シミュレータの使用方法)
6. [トラブルシューティング](#トラブルシューティング)

---

## クイックスタート
バッチファイルを走らせると、必要な環境の有無の確認と必要な環境のインストールを行うように作ってあります。

> [!NOTE]
> **Windows専用**: このセットアップスクリプトは Windows 環境でのみ動作します。
> 
> Linux / macOS ニキは [手動セットアップ](#手動セットアップ) を参照してください。

リポジトリをクローンした後、以下のコマンドを実行するだけで環境構築とビルドが完了します：

### 完全自動セットアップ（推奨）

```batch
git clone <repository-url>
cd research-tbp
tools\setup.bat --unattended --with-oneapi
```

`--unattended` オプションにより、すべての依存関係が自動的にインストールされます。
`--with-oneapi` オプションによりIntel oneAPI のインストールが有効化されます（インストールに10-30分かかりますよドンマイ）。

### 完全自動セットアップ（gcc代替）

```batch
tools\setup.bat --unattended --with-oneapi
```

### 対話モード（何がインストールされるか確認しながらやりたいとき）

```batch
git clone <repository-url>
cd research-tbp
tools\setup.bat
```

### セットアップオプション

| オプション | 説明 |
|-----------|------|
| `--unattended`, `-y` | 非対話モード（すべて自動インストール） |
| `--with-oneapi` | Intel oneAPI をインストール（推奨、10-30分） |
| `--skip-python` | Python環境セットアップをスキップ |
| `--skip-test` | ビルド後のテスト実行をスキップ |
| `--help`, `-h` | ヘルプを表示 |

### セットアップスクリプトが行うこと

1. 依存関係のチェック（Git, CMake, Ninja, コンパイラ）
2. 不足しているツールの自動インストール（winget経由）
3. vcpkg と Boost のインストール（`external/` ディレクトリに）
4. `BOOST_ROOT` 環境変数の永続的な登録
5. 最適なコンパイラの自動選択とプロジェクトのビルド
6. Python依存パッケージのインストール（オプション）
7. テスト実行によるビルド検証（オプション）

---

## 必要な環境

### 必須

| ツール | バージョン | インストール方法 |
|--------|-----------|-----------------|
| Git | 任意 | https://git-scm.com/ |
| CMake | 3.15以上 | https://cmake.org/ |

### コンパイラ（いずれか1つ）

| コンパイラ | 推奨度 | インストール方法 |
|-----------|-------|-----------------|
| Intel oneAPI + Ninja | 推奨 | [Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html) |
| MinGW-gcc | 代替 | `winget install -e --id mingw.mingw-w64-ucrt-x86_64` |

### オプション

| ツール | 用途 | インストール方法 |
|--------|------|-----------------|
| Ninja | 高速ビルド | `winget install Ninja-build.Ninja` |

---

## 手動セットアップ

自動セットアップを使用しない場合は、以下の手順で環境を構築できます：

### 1. Boost のインストール

**方法A: vcpkg を使用（推奨）**
```batch
cd external
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg
bootstrap-vcpkg.bat -disableMetrics
vcpkg install boost:x64-windows
```

**方法B: 手動インストール**
1. [Boost公式サイト](https://www.boost.org/)からダウンロード
2. 任意の場所に展開

### 2. 環境変数の設定

`BOOST_ROOT` をユーザー環境変数に設定します：

```batch
rem vcpkg を使用した場合
setx BOOST_ROOT "C:\path\to\TBP\external\vcpkg\installed\x64-windows"

rem 手動インストールの場合
setx BOOST_ROOT "C:\path\to\boost_1_xx_0"
```

### 3. ビルド

**Intel oneAPI + Ninja の場合:**
```batch
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
cmake --preset=ninja-intel
cmake --build build
```

**MinGW-gcc の場合:**
```batch
cmake --preset=mingw-gcc
cmake --build build
```

---

## ディレクトリ構成

```
TBP/
├── CMakeLists.txt            # ルートのCMake設定
├── CMakePresets.json         # CMakeビルドプリセット
├── README.md                 # このファイル
│
├── scripts/                  # セットアップスクリプト
│   ├── setup.bat             # メインセットアップスクリプト
│   ├── check_dependencies.bat # 依存関係チェック
│   └── install_vcpkg_boost.bat # vcpkg/Boostインストール
│
├── external/                 # 外部依存ライブラリ
│   └── vcpkg/                # vcpkg (setup.bat実行後に作成)
│
├── apps/                     # シミュレーションアプリケーション
│   ├── SALI2d/               # 2D SALI解析
│   ├── SALI3d/               # 3D SALI解析
│   ├── PoincareMap/          # ポアンカレ断面
│   └── ...
│
├── lib/                      # 共有ライブラリ
│   └── include/              # ヘッダファイル
│       ├── rtbp.hpp          # CR3BP計算クラス
│       ├── crtbp.hpp         # CR3BP追加機能
│       └── utils.hpp         # ユーティリティ関数
│
├── configs/                  # シミュレーション設定ファイル
├── data/                     # 出力データ
├── test/                     # テストコード
└── build/                    # ビルド出力 (gitignore)
```

---

## シミュレータの使用方法

ビルド後、実行ファイルは `build/bin/` に生成されます：

```batch
rem 例: SALI2dの実行
build\bin\Debug\SALI2d.exe configs\sali2d\config.toml

rem 例: ポアンカレマップの生成
build\bin\Debug\PoincareMap.exe configs\poincare\config.toml
```

---

## トラブルシューティング

### BOOST_ROOT が見つからない

```
CMake Error: Could NOT find Boost
```

**解決策:**
1. 新しいターミナルを開き、`echo %BOOST_ROOT%` で環境変数を確認
2. 環境変数が設定されていない場合は `scripts\setup.bat` を再実行

### Ninja が見つからない

```
CMake Error: CMake was unable to find a build program corresponding to "Ninja"
```

**解決策:**
1. `winget install Ninja-build.Ninja` でインストール
2. または `mingw-gcc-make` プリセットを使用: `cmake --preset=mingw-gcc-make`

### oneAPI 環境が設定されていない

```
'icx-cl' is not recognized as an internal or external command
```

**解決策:**
1. `call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"` を実行
2. または MinGW-gcc プリセットを使用

---

## 参考リンク

- [Intel oneAPI インストールガイド](https://dexall.co.jp/articles/?p=1968)
- [CMake + Intel コンパイラ設定](https://qiita.com/Sego-don/items/f5c8adf853c4badf8171)

