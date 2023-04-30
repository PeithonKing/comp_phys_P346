# Computational Physics Lab (P346) Assignments

(readme last updated 30th November 2022)

## Assignments

1. [**Assignment 1:** General Programming](1.assign1_gen.ipynb)
2. [**Assignment 2:** Random Number Generation](2.assign2_rand.ipynb)
3. [**Assignment 3:** Solution of Linear Equations](3.assign3_lineq.ipynb)
4. [**Assignment 4:** Roots of nonlinear equations, polynomials and data-fitting](4.assign4_rootfit.ipynb)
5. [**Assignment 5:** Numerical Integration](5.assign5_integ.ipynb)
6. [**Assignment 6:** Differential Equations and eigen values](6.assign6_difqev.ipynb)

## Exams and Projects

- [**Midsem Exam**](midsem.ipynb)
- [**DIY Project**](DIY%20Project/main.pdf)
- [**Endsem Exam**](endsem.ipynb)

## Lab Works

1. [**11th August:** Random Number Generation](lab_11_08.ipynb)
2. [**17th August:** Gauss-Jordan Elimination](lab_17_08.ipynb)
3. [**24th August:** LU Doolittle's and Cholesky Decomposition](lab_24_08.ipynb)
4. [**31st August:** Gauss-Seidel Method](lab_31_08.ipynb)
5. [**7th September:** Bisection, Regula-Falsi and Newton-Raphson Method](lab_07_09.ipynb)
6. [**14th September:** Laguerre method and deflation](lab_14_09.ipynb)
7. [**15th September:** Interpolation](lab_15_09.ipynb)
8. [**12th October:** Data fitting](lab_12_10.ipynb)
9. [**13th October:** Numerical Integration](lab_13_10.ipynb)
10. [**26th October:** Monte Carlo Integration](lab_26_10.ipynb)
11. [**27th October:** Runge-Kutta Methods](lab_27_10.ipynb)
12. [**9th November:** Partial Differential Equations](lab_9_11.ipynb)
13. [**10th November:** Eigen Values and Eigen Vectors](lab_10_11.ipynb)
14. [**16th November:** Radioactivity](lab_16_11.ipynb)

## Deployment (How to run codes in this repository)

&emsp; Hope you have Python and Git Installed. I have used Python 3.9.13 version which is easily available in the windows store. If not, you can install it from [this](https://www.python.org/downloads/release/python-390/) link.

1. **Step 1:** Clone the repository

    ```powershell
    git clone https://github.com/PeithonKing/comp_phys_P346.git
    ```

2. **Step 2:** Go into the directory

    ```powershell
    cd comp_phys_P346
    ```

3. **Step 3:** Make a virtual environment

    ```powershell
    python3 -m vitualenv venv
    ```

    if you don't have virtualenv installed, you can install it by the following command before doing the above command

    ```powershell
    pip install virtualenv
    ```

4. **Step 4:** Activate the virtual environment
    - in Windows:

    ```powershell
    venv\Scripts\activate
    ```

    - in Linux / MacOS:

    ```powershell
    source venv/bin/activate
    ```

5. **Step 5:** Install all the packages

    ```powershell
    pip install -r requirements.txt
    ```

6. **Step 6:** Open the notebook

    ```powershell
    jupyter notebook
    ```

    or if jupyter is not added to the system path, you can use

    ```powershell
    python -m notebook
    ```

7. **Step 7:** Run the notebooks
