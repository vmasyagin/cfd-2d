# cfd-2d
Automatically exported from code.google.com/p/cfd-2d

# Основная цель проекта

Создание двумерного и трехмерного CFD-кодов для моделирования турбулентных течений сжимаемого газа на основе схем высокого порядка точности. В качестве расчетной методики будут использованы явные и неявные конечно-объемные методики. Также предполагается создание кодов на основе разрывного метода Галеркина с использованием явных и неявных методов дискретизации уравнений Навье-Стокса по времени.


# Текущие задачи

Обобщение и реализация на многомерные уравнения газовой динамики энтропийного ограничителя наклона (http://www.mathnet.ru/links/086f7241e147b53c9302f01c250c91ab/ipmp2689.pdf). 

Реализация трехмерной неявной методики на основе метода Галеркина с разрывными базисными функциями с учетом эффектов вязкости.

Реализация лимитера Кокбурна для трехмерных уравнений газовой динамики.

Реализация лимитера для разрывного метода Галеркина на основе усреднения решения для трехмерных уравнений газовой динамики (http://www.mathnet.ru/links/98d3b3ee1ef7a4653baa8f55eb63a3b1/mm3970.pdf).

Верификация программного кода.

Реализация поддержки форматов хранения данных популярных CFD-решателей.

Реализация RANS-моделей турбулентности

Участники проекта
Сотрудники и аспиранты кафедры прикладной математики, дифференциальных уравнений и теоретической механики Мордовского государственного университета им. Н.П.Огарева.
