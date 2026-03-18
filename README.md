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

### ビルドシステム

| ツール | 用途 | インストール方法 |
|--------|------|-----------------| 
| Ninja | 高速ビルド（推奨） | `winget install Ninja-build.Ninja` |
| Visual Studio | Intel oneAPI使用時に必要 | [Visual Studio](https://visualstudio.microsoft.com/) |

---

### C++ 外部ライブラリ

| ライブラリ | 用途 | 導入方法 |
|-----------|------|---------| 
| Boost (≥1.71) | ODE数値積分 | vcpkg または手動（`setup.bat` で自動導入） |
| toml++ (v3.4.0) | TOML設定ファイルパーサー | CMake FetchContent（自動取得） |
| Eigen3 (3.4.0) | 線形代数・固有値計算 | CMake FetchContent（自動取得） |
| OpenMP | 並列計算 | コンパイラに付属 |

## 手動セットアップ

自動セットアップを使用しない場合は、以下の手順で環境を構築：

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

`BOOST_ROOT` をユーザー環境変数に設定：

```batch
rem vcpkg を使用した場合
setx BOOST_ROOT "C:\path\to\TBP\external\vcpkg\installed\x64-windows"

rem 手動インストールの場合
setx BOOST_ROOT "C:\path\to\boost_1_xx_0"
```

### 3. ビルド

**Intel oneAPI + Visual Studio の場合:**
```batch
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
cmake --preset=vs-intel
cmake --build build
```

**Intel oneAPI + Ninja の場合:**
```batch
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
cmake --preset=ninja-intel
cmake --build build
```

**MinGW-gcc + Ninja の場合:**
```batch
cmake --preset=mingw-gcc
cmake --build build
```

**MinGW-gcc + MinGW Makefiles の場合（Ninja不要）:**
```batch
cmake --preset=mingw-gcc-make
cmake --build build
```

---

## ディレクトリ構成

```
TBP/
├── CMakeLists.txt            # ルートのCMake設定
├── CMakePresets.json         # CMakeビルドプリセット
├── README.md                 # このファイル
├── build_all.bat             # 全アプリ一括ビルドスクリプト
├── test.bat                  # テスト実行スクリプト
│
├── tools/                    # セットアップ・解析ツール
│   ├── setup.bat             # メインセットアップスクリプト
│   ├── check_dependencies.bat # 依存関係チェック
│   ├── install_vcpkg_boost.bat # vcpkg/Boostインストール
│   └── *.py                  # マニフォールド解析等のPythonツール
│
├── scripts/                  # Python解析スクリプト群
│   ├── saliplotter.py        # SALI可視化ツール
│   ├── zvc_2d_interactive.py # ZVC 2Dインタラクティブ描画
│   ├── zvc_3d_interactive.py # ZVC 3Dインタラクティブ描画
│   ├── crtbp.py              # CR3BP Python実装
│   ├── poincare_map.py       # ポアンカレ断面描画
│   ├── requirements.txt      # Python依存パッケージ
│   └── ...
│
├── external/                 # 外部依存ライブラリ
│   └── vcpkg/                # vcpkg (setup.bat実行後に作成)
│
├── apps/                     # シミュレーションアプリケーション
│   ├── SALI2d/               # 2D SALI解析
│   ├── SALI3dV3/             # 3D SALI解析 v3
│   ├── SaliOrbPlane/         # SALI軌道面解析
│   ├── PoincareMap/          # ポアンカレ断面
│   ├── TrajectoryCalc/       # 軌道計算
│   ├── TrajectorySali/       # 軌道SALI計算
│   ├── PeriodicOrbitAnalysis/# 周期軌道解析
│   ├── StableManifoldComputation/ # 安定マニフォールド計算
│   ├── OrbitManifoldByJacobi/# ヤコビ積分によるマニフォールド
│   ├── ImpulseOrbitSearch/   # インパルス軌道探索
│   ├── AsteroidOrbitIntersection/ # 小惑星軌道交差計算
│   ├── JacobiIntegralCalc/   # ヤコビ積分計算
│   ├── JacobiVelocityPlot/   # ヤコビ速度プロット
│   ├── ChaosIndicator/       # カオス指標計算
│   ├── CalcMethodCompare/    # 計算手法比較
│   ├── CoordTransformVerify/ # 座標変換検証
│   └── sticky_orbit_search/  # 粘着軌道探索
│
├── lib/                      # 共有ライブラリ
│   ├── CMakeLists.txt
│   ├── include/              # ヘッダファイル
│   │   ├── rtbp.hpp          # CR3BP計算コアクラス
│   │   ├── crtbp.hpp         # CR3BP追加機能
│   │   ├── utils.hpp         # ユーティリティ関数
│   │   ├── vector3d.hpp      # 3Dベクトルクラス
│   │   ├── periodic_orbit.hpp # 周期軌道計算
│   │   ├── periodic_orbit_impl.hpp # 周期軌道実装
│   │   ├── periodic_orbit_stability.hpp # 周期軌道安定性解析
│   │   ├── periodic_orbit_manifold.hpp  # 周期軌道マニフォールド
│   │   ├── lyapunov_initial.hpp # Lyapunov軌道初期値
│   │   └── continuation.hpp  # 連続法
│   └── src/                  # ソースファイル
│
├── configs/                  # シミュレーション設定ファイル (TOML)
├── configdata/               # 設定データテンプレート
├── docs/                     # ドキュメント
├── data/                     # 出力データ
├── test/                     # テストコード
├── __past/                   # 過去のコード・参考資料
└── build/                    # ビルド出力 (gitignore)
```

---

## シミュレータの使用方法

ビルド後、実行ファイルは `build/bin/` に生成されます：

```batch
rem 例: SALI2dの実行
build\bin\SALI2d.exe configs\3D_crtbp_SALI\config.toml

rem 例: ポアンカレマップの生成
build\bin\PoincareMap.exe configs\poincare_map\config.toml
```

### 全アプリの一括ビルド

```batch
build_all.bat
```

### Python解析スクリプトの使用

```batch
cd scripts
pip install -r requirements.txt
python saliplotter.py
```

---

## トラブルシューティング

### BOOST_ROOT が見つからない

```
CMake Error: Could NOT find Boost
```

**解決策:**
1. 新しいターミナルを開き、`echo %BOOST_ROOT%` で環境変数を確認
2. 環境変数が設定されていない場合は `tools\setup.bat` を再実行

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
