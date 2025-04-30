#ifndef GNUPLOT_HPP
#define GNUPLOT_HPP

#include <memory>
#include <mutex>
#include <queue>
#include <string>

class GnuplotCommand {
   public:
    // コンストラクタ
    explicit GnuplotCommand(const std::string& gnuplotPath = "gnuplot");

    // デストラクタ
    ~GnuplotCommand();

    // コピーコンストラクタとコピー代入演算子を削除
    GnuplotCommand(const GnuplotCommand&) = delete;
    GnuplotCommand& operator=(const GnuplotCommand&) = delete;

    // ムーブコンストラクタとムーブ代入演算子
    GnuplotCommand(GnuplotCommand&&) noexcept;
    GnuplotCommand& operator=(GnuplotCommand&&) noexcept;

    // コマンド実行メソッド
    template <typename... Args>
    void command(const char* format, Args&&... args);

    // フラッシュメソッド
    void flush();

   private:
    std::string gnuplotPath_;
    std::unique_ptr<FILE, decltype(&pclose)> fp_;
    mutable std::mutex mutex_;
    std::queue<std::string> commandQueue_;

    void openGnuplotProcess();
    void closeGnuplotProcess();

    template <typename... Args>
    std::string formatString(const char* format, Args&&... args);
};

// テンプレートメソッドの宣言
// template<typename... Args>
// void GnuplotCommand::command(const char* format, Args&&... args);

// template<typename... Args>
// std::string GnuplotCommand::formatString(const char* format, Args&&... args);

#endif  // GNUPLOT_HPP