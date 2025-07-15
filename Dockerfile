# Dockerfile using renv

# 1. 选择一个基础镜像
# 推荐使用 rocker 项目提供的镜像
FROM rocker/r-ver:4.4.1

# 2. 安装系统依赖
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
    # --- 基础依赖 ---
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    # --- 根据R包编译错误逐步添加 ---
    cmake \
    libcairo2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    pandoc \
    perl \
    libharfbuzz-dev \
    libfribidi-dev \
    libtiff5-dev \
    libjpeg-dev \
    # --- 清理apt缓存以减小镜像体积 ---
    && rm -rf /var/lib/apt/lists/*

# 3. 安装 renv 包
RUN R -e "install.packages('renv')"

# 4. 拷贝 renv 相关文件到工作目录
WORKDIR /app
COPY renv.lock .
COPY .Rprofile .
COPY renv/activate.R renv/

# 本地下载gPCA包
# 解决gPCA包安装失败
COPY ./utils/gPCA_1.0.tar.gz .
RUN R -e "install.packages('gPCA_1.0.tar.gz', repos = NULL, type = 'source')"

# 强制C++编译器使用C++14标准
# 解决fgsea包编译失败
ENV CXX14FLAGS="-O2 -Wall -pedantic -std=c++14"

# 5. 利用 renv::restore() 恢复环境
# 利用Docker的缓存机制，只要renv.lock不变，此步骤不会重复执行
RUN R -e "renv::restore()"

# 6. 将本地的其他文件（R脚本、数据等）复制到容器的工作目录中
COPY . .

# 7. 定义容器启动时执行的命令
# 将 Rscript 设置为入口点程序
ENTRYPOINT ["Rscript"]
# 将 my_script.R 设置为默认要运行的脚本
CMD ["test.R"]

# 8. 限制最大并行数
ENV OPENBLAS_NUM_THREADS=8
