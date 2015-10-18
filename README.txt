README

Системные требования:
* python 2.7.*
* python-bio

1) Собрать файлы, написанные на Cython. Из корневой директории выполнить:
python setup.py build_ext --inplace

2) Запустить код программы с профайлером:
python -m cProfile -o timeStats.profile Main.py

3) Выбрать режим работы.

4) После работы программы можно проанализировать статистику.
python -m pstats timeStats.profile

5) Выбрать функцию, для которой просматривается статистика, пример:
stats test_debruijn

6) Для выхода ввести:
quit
