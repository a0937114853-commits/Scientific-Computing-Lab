# Scientific Computing Lab - Poisson Equation Solvers

本專案為科學計算（Scientific Computing）相關作業與程式碼實作，主要使用 **MATLAB** 實作二維帕松方程式（Poisson Equation）的數值求解器，涵蓋 5 點差分法（5-point Laplacian）、矩形網格非均勻間距推廣、整數網格自動修正機制，以及達到 4 次方精確度的 9 點高階差分法（9-point Laplacian）。

相關的數學推導與數值分析詳見 [Scientific Computing Hw3.pdf](Scientific%20Computing%20Hw3.pdf)。

---

## 專案目錄與檔案說明

專案中的 MATLAB 程式碼主要對應以下數值方法與實驗：

* **`poisson_rect_general.m`**：
  * 實作二維帕松方程式在一般矩形域 $\Omega = [a_x, b_x] \times [a_y, b_y]$ 上的 5 點差分法。
  * 支援非均勻網格間距（$\Delta x \neq \Delta y$），透過 Kronecker 產品有效建構大型稀疏矩陣 $A$。
* **`poisson_integer_grid_fix.m`**：
  * 針對非正方形矩形域給定目標間距 $h_{target}$ 時，可能導致網格點數非整數的「奇特現象（Strange Phenomenon）」進行修正。
  * 透過 `round()` 機制自動調整網格數，確保邊界切合並維持 $O(h^2)$ 精度。
* **`convergence_study_square.m`**：
  * 在正方形區域上執行網格加密研究（Grid Refinement Study），依序測試 $m = 20, 40, 80, 160$。
  * 計算最大誤差（$L_\infty$ norm）並透過 Log-Log 圖驗證其理論上的 **第二階收斂精度 $O(h^2)$**。
* **`poisson_9point_solver.m`**：
  * 實作高效的 **9 點拉普拉斯算子（9-point Laplacian）**，並對源項進行 $f + \frac{h^2}{12}\nabla^2 f$ 的修正。
  * 透過邊界補償項與 Kronecker 矩陣組合，成功在數值實驗中驗證並達到 **第四階收斂精度 $O(h^4)$**。

---

## 數學模型與測試問題

考慮二維帕松方程式與 Dirichlet 邊界條件：
$$\nabla^2 u(x, y) = \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} = f(x, y), \quad \text{in } \Omega$$

測試問題中所使用的解析解與源項（Source Term）設定為：
$$f(x, y) = 1.25 \exp(x + y/2)$$
$$u_{true}(x, y) = \exp(x + y/2)$$

---

## 執行方式

1. 確認電腦已安裝 **MATLAB**。
2. 將此專案複製（Clone）或下載至本地端。
3. 在 MATLAB 中打開對應的 `.m` 檔（例如 `convergence_study_square.m` 或 `poisson_9point_solver.m`），點擊 **Run** 即可直接執行並查看結果與誤差收斂圖表。

---
*Author: Chao-Yang Zhan*