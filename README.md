# TBP - Three Body Problem Simulation

円制限三体問題（CR3BP）の数値計算・解析リポジトリ。

## 目次

1. 最短セットアップ（初回はこれだけ）
2. セットアップの詳細
3. 手動セットアップ
4. 実行方法
5. テストと再ビルド
6. リポジトリ構成
7. トラブルシュート

## 1. 最短セットアップ（初回はこれだけ）

この節の手順だけで、何も環境がない Windows マシンからビルドまで進められる設計になっている。

### 前提

- OS: Windows 10/11
- PowerShell か cmd を使用
- インターネット接続あり

### 手順

```batch
git clone <repository-url>
cd TBP
tools\setup.bat --unattended --python-venv
```

`setup.bat` が以下を自動で処理する。

- 依存ツール確認（Git, CMake, Ninja, oneAPI, MinGW, Windows SDK）
- 不足依存の自動導入（winget が正常な場合）
- `external/vcpkg` と Boost の導入
- CMake configure/build
- Python 依存の `.venv` への導入（`--python-venv` 指定時）

このコマンドは oneAPI をデフォルトで導入対象にする。

### oneAPI を導入しない場合

```batch
tools\setup.bat --unattended --without-oneapi --python-venv
```

### 最短確認

```batch
dir build\bin
```

`build/bin` に実行ファイルが生成されていればセットアップ完了。

## 2. セットアップの詳細

### setup.bat の主なオプション

| オプション | 内容 |
|---|---|
| `--unattended`, `-y` | 非対話モード。確認プロンプトなしで自動処理する |
| `--with-oneapi` | Intel oneAPI の導入を有効化する（デフォルト） |
| `--without-oneapi` | Intel oneAPI の導入を無効化する |
| `--skip-python` | Python 依存の導入をスキップする |
| `--python-venv` | Python 依存を `.venv` に導入する |
| `--skip-test` | ビルド後のテスト実行をスキップする |
| `--help`, `-h` | ヘルプを表示する |

### 推奨実行パターン

### 完全自動（推奨）

```batch
tools\setup.bat --unattended --python-venv
```

### oneAPI を無効化する場合

```batch
tools\setup.bat --unattended --without-oneapi --python-venv
```

### 対話モード

```batch
tools\setup.bat
```

### コンパイラ選択ロジック

`setup.bat` は利用可能な環境に応じて自動選択する。

- oneAPI + Ninja が使える: `ninja-intel`
- oneAPI 不可かつ MinGW + Ninja が使える: `mingw-gcc`
- Ninja なしだが MinGW が使える: `mingw-gcc-make`

### winget について

自動インストールは `winget` が正常動作することを前提にする。

- `check_dependencies.bat` は `winget --info` で健全性を確認する
- `winget` が壊れている場合、App Installer の修復を試行する
- 修復できない場合、自動インストールを無効化して警告を表示する

## 3. 手動セットアップ

ネットワーク制限やポリシーで自動導入できない場合の手順。

### 1) 必須ツール導入

- Git
- CMake 3.15+
- コンパイラ（以下のどちらか）
- Intel oneAPI + Ninja
- MinGW-gcc

### 2) Boost 導入（vcpkg）

```batch
cd external
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg
bootstrap-vcpkg.bat -disableMetrics
vcpkg install boost-odeint:x64-windows
```

### 3) BOOST_ROOT 設定

```batch
setx BOOST_ROOT "C:\path\to\TBP\external\vcpkg\installed\x64-windows"
```

### 4) CMake configure/build

#### Intel oneAPI + Ninja

```batch
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
cmake --preset=ninja-intel
cmake --build build
```

#### MinGW + Ninja

```batch
cmake --preset=mingw-gcc
cmake --build build
```

#### MinGW + Makefiles（Ninjaなし）

```batch
cmake --preset=mingw-gcc-make
cmake --build build
```

### 5) Python 依存導入（任意）

```batch
python -m venv .venv
.venv\Scripts\python -m pip install -r scripts\requirements.txt
```

## 4. 実行方法

ビルド後の実行ファイルは `build/bin` に生成される。

### 例: SALI2d

```batch
build\bin\SALI2d.exe configs\3D_crtbp_SALI\3DSALIconfig_sample.toml
```

### 例: PoincareMap

```batch
build\bin\PoincareMap.exe configs\poincare_map\poincare_map_config_sample.toml
```

## 5. テストと再ビルド

### 再ビルド

```batch
cmake --build build
```

### テスト実行

```batch
ctest --test-dir build --output-on-failure -C Debug
```

### 一括ビルド

```batch
build_all.bat
```

## 6. リポジトリ構成

主要ディレクトリのみ記載。

```text
TBP/
├── apps/                 # 各シミュレーション実行アプリ
├── lib/                  # 共通ライブラリ
├── configs/              # 各アプリ用 TOML 設定
├── scripts/              # Python 解析スクリプト
├── tools/                # セットアップ/補助バッチ
├── test/                 # C++ テスト
├── external/             # vcpkg など外部依存
├── CMakeLists.txt
├── CMakePresets.json
└── README.md
```

## 7. トラブルシュート

### `winget` が使えない

症状例:

- `winget --info` が失敗する
- `The file cannot be accessed by the system`

対応:

1. Microsoft Store から App Installer を更新/再インストール
2. 新しいターミナルを開いて `winget --info` を再確認
3. それでも不可なら手動セットアップに切り替える

### Boost が見つからない

症状例:

```text
CMake Error: Could NOT find Boost
```

対応:

1. `echo %BOOST_ROOT%` で環境変数確認
2. `external\vcpkg\installed\x64-windows\include\boost` の存在確認
3. `tools\setup.bat` を再実行

### Ninja が見つからない

症状例:

```text
CMake Error: CMake was unable to find a build program corresponding to "Ninja"
```

対応:

1. `winget install Ninja-build.Ninja`
2. もしくは `cmake --preset=mingw-gcc-make` を使用

### oneAPI が見つからない

症状例:

```text
'icx-cl' is not recognized as an internal or external command
```

対応:

1. `call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"`
2. もしくは MinGW プリセットでビルド
