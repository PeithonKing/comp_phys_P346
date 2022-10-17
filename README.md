# Computational Physics Lab (P346) Assignments

(readme last updated on 12th October 2022)

## Assignments

1. [**Assignment 1:** General Programming](1.assign1_gen.ipynb)
2. [**Assignment 2:** Random Number Generation](2.assign2_rand.ipynb)
3. [**Assignment 3:** Solution of Linear Equations](3.assign3_lineq.ipynb)
4. 

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
10. 

## Deployment (How to run codes in this repository)

&emsp; Hope you have Python and Git Installed. I have used Python 3.9.13 version which is easily available in the windows store. If not, you can install it from [this](https://www.python.org/downloads/release/python-390/) link.

1. **Step 1:** Clone the repository

    ```bash
    git clone https://github.com/PeithonKing/comp_phys_P346.git
    ```

2. **Step 2:** Enter the directory
    ```bash
    cd comp_phys_P346
    ```

3. **Step 3:** Make a virtual environment
    ```bash
    python3 -m vitualenv venv
    ```

    if you don't have virtualenv installed, you can install it by
    ```bash
    pip install virtualenv
    ```

4. **Step 4:** Activate the virtual environment
    * in Windows:
    ```bash
    venv\Scripts\activate
    ```
    * in Linux / MacOS:
    ```bash
    source venv/bin/activate
    ```

5. **Step 5:** Install all the packages

    ```bash
    pip install -r requirements.txt
    ```

6. **Step 6:** Open the notebook
    ```bash
    jupyter notebook
    ```
    or if jupyter is not added to the system path, you can use
    ```bash
    python -m notebook
    ```

7. **Step 7:** Run the notebooks
