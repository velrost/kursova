import math
from numpy import *
from sympy import *
from math import *



def cur_func(x1, x2):
    global func_usage
    func_usage = func_usage + 1
    return Stepen_func(x1, x2)

def Stepen_func(x1, x2):
    return (10 * (x1 - x2) ** 2 + (x1 - 1) ** 2) ** 1/4


def left_derivative(point, h, var_to_count, cur_func_value):
    if var_to_count == "x1":
        return (cur_func_value - cur_func(point[0] - h, point[1])) / h
    if var_to_count == "x2":
        return (cur_func_value - cur_func(point[0], point[1] - h)) / h


def right_derivative(point, h, var_to_count, cur_func_value):
    if var_to_count == "x1":
        return (cur_func(point[0] + h, point[1]) - cur_func_value) / h
    if var_to_count == "x2":
        return (cur_func(point[0], point[1] + h) - cur_func_value) / h


def central_derivative(point, h, var_to_count):
    if var_to_count == "x1":
        return (cur_func(point[0] + h, point[1]) - cur_func(point[0] - h, point[1])) / (
            2 * h
        )
    if var_to_count == "x2":
        return (cur_func(point[0], point[1] + h) - cur_func(point[0], point[1] - h)) / (
            2 * h
        )


def find_gradient(point, cur_func_value):

    return [
        left_derivative(point, h, "x1", cur_func_value),
        left_derivative(point, h, "x2", cur_func_value),
    ]


    #return [right_derivative(point, h,"x1", cur_func_value), right_derivative(point, h,"x2", cur_func_value)]
    #return [central_derivative(point, h,"x1"), central_derivative(point, h,"x2")]


def find_norm(thing):
    return abs(math.sqrt(thing[0] * thing[0] + thing[1] * thing[1]))


def find_delta_lambda(norm_point, norm_step):
    global parametr_Svena
    try:
        step = parametr_Svena * (norm_point / norm_step)
    except ZeroDivisionError:
        step = parametr_Svena
    return step


def System_End_Criteria(cur_point, cur_func_value, prevoius_point, prevoius_func_value):
    point_difference = [
        cur_point[0] - prevoius_point[0],
        cur_point[1] - prevoius_point[1],
    ]
    point_error = find_norm(point_difference) / find_norm(prevoius_point)

    func_error = abs(cur_func_value - prevoius_func_value)
    return [point_error, func_error]


def Sven_algorithm(delta_lambda, point, func_value, gradient):
    local_delta_lambda = 0
    cur_lambda = 0
    plus_lambda = delta_lambda
    minus_lambda = -delta_lambda
    plus_point = [
        point[0] - plus_lambda * gradient[0],
        point[1] - plus_lambda * gradient[1],
    ]
    minus_point = [
        point[0] - minus_lambda * gradient[0],
        point[1] - minus_lambda * gradient[1],
    ]

    print(
        f"\n\n АЛГОРИТМ СВЕНА \n\n Start point: {point},  {func_value}\nDelta Lambda:  {delta_lambda} \nAnti-Gradient: {[-gradient[0],-gradient[1]]} \n"
    )

    plus_func = cur_func(plus_point[0], plus_point[1])
    minus_func = cur_func(minus_point[0], minus_point[1])

    print(
        f"Plus delta: {plus_point},  {plus_func}\nMinus Delta {minus_point},  {minus_func}\n"
    )

    global next_func_value
    global cur_func_value
    cur_func_value = func_value

    print(f"Начальное значение функции: {cur_func_value}")
    if plus_func > func_value and minus_func > func_value:
        while plus_func > func_value and minus_func > func_value:
            plus_lambda = plus_lambda / 2
            minus_lambda = minus_lambda / 2
            plus_point = [
                point[0] - plus_lambda * gradient[0],
                point[1] - plus_lambda * gradient[1],
            ]
            minus_point = [
                point[0] - minus_lambda * gradient[0],
                point[1] - minus_lambda * gradient[1],
            ]

            plus_func = cur_func(plus_point[0], plus_point[1])
            minus_func = cur_func(minus_point[0], minus_point[1])
            print(
                f"NEW plus delta: {plus_point},  {plus_func}\nNEW Minus Delta {minus_point},  {minus_func}\n"
            )
    if plus_func < func_value:
        local_delta_lambda = plus_lambda
        point = plus_point
        next_func_value = plus_func
        func_value = plus_func
    elif minus_func < func_value:
        local_delta_lambda = minus_lambda
        point = minus_point
        next_func_value = minus_func

    print(f"Теперешнее значение функции {next_func_value}\n")

    points_arr = []
    func_value_arr = []
    lambda_valu_arr = []
    points_arr.append(point)
    func_value_arr.append(next_func_value)
    lambda_valu_arr.append(local_delta_lambda)

    k = 1
    while next_func_value < cur_func_value:
        cur_lambda = cur_lambda + local_delta_lambda * 2 ** k
        new_point = [
            point[0] - cur_lambda * gradient[0],
            point[1] - cur_lambda * gradient[1],
        ]
        points_arr.append(new_point)

        lambda_valu_arr.append(cur_lambda)
        cur_func_value = next_func_value
        next_func_value = cur_func(new_point[0], new_point[1])
        func_value_arr.append(next_func_value)

        print(f"Новая точка: {new_point}")
        print(f"Теперешнее значение функции {next_func_value}\n")
        k = k + 1

    cur_lambda = cur_lambda - local_delta_lambda * k
    new_point = [
        point[0] - cur_lambda * gradient[0],
        point[1] - cur_lambda * gradient[1],
    ]
    points_arr.append(new_point)
    lambda_valu_arr.append(cur_lambda)
    print(f"Новая точка: {new_point}")
    next_func_value = cur_func(new_point[0], new_point[1])
    func_value_arr.append(next_func_value)
    print(f"Теперешнее значение функции {next_func_value}\n")

    return [
        [
            points_arr[len(points_arr) - 4],
            func_value_arr[len(func_value_arr) - 4],
            lambda_valu_arr[len(points_arr) - 4],
        ],
        [
            points_arr[len(points_arr) - 3],
            func_value_arr[len(func_value_arr) - 3],
            lambda_valu_arr[len(points_arr) - 3],
        ],
        [
            points_arr[len(points_arr) - 1],
            func_value_arr[len(points_arr) - 1],
            lambda_valu_arr[len(points_arr) - 1],
        ],
    ]


def DSK_Pauela(left_border_arr, middle_arr, right_border_arr, gradient):
    print("\nДСК-ПАУЭЛА\n")
    delta_x = [
        middle_arr[0][0] - left_border_arr[0][0],
        middle_arr[0][1] - left_border_arr[0][1],
    ]
    try:
        star_x = [
            middle_arr[0][0]
            + (
                (delta_x[0] * (left_border_arr[1] - right_border_arr[1]))
                / (2 * (left_border_arr[1] - 2 * middle_arr[1] + right_border_arr[1]))
            ),
            middle_arr[0][1]
            + (
                (delta_x[1] * (left_border_arr[1] - right_border_arr[1]))
                / (2 * (left_border_arr[1] - 2 * middle_arr[1] + right_border_arr[1]))
            ),
        ]
    except ZeroDivisionError:
        star_x = [middle_arr[0], middle_arr[1]]
        return star_x

    star_func_value = cur_func(star_x[0], star_x[1])
    print(f"Point : {star_x}\n Function: {star_func_value}")

    global tochnost_odn
    while abs(left_border_arr[1] - right_border_arr[1]) > tochnost_odn:
        if left_border_arr[1] < right_border_arr[1]:
            right_border_arr = middle_arr
        else:
            left_border_arr = middle_arr

        middle_arr = [star_x, star_func_value]
        a_1 = [
            (middle_arr[1] - left_border_arr[1])
            / (middle_arr[0][0] - left_border_arr[0][0]),
            (middle_arr[1] - left_border_arr[1])
            / (middle_arr[0][1] - left_border_arr[0][1]),
        ]
        a_2 = [
            (1 / (right_border_arr[0][0] - middle_arr[0][0]))
            * (
                (
                    (right_border_arr[1] - left_border_arr[1])
                    / (right_border_arr[0][0] - left_border_arr[0][0])
                )
                - a_1[0]
            ),
            (1 / (right_border_arr[0][1] - middle_arr[0][1]))
            * (
                (
                    (right_border_arr[1] - left_border_arr[1])
                    / (right_border_arr[0][1] - left_border_arr[0][1])
                )
                - a_1[1]
            ),
        ]
        star_x = [
            ((left_border_arr[0][0] + middle_arr[0][0]) / 2) - (a_1[0] / (2 * a_2[0])),
            ((left_border_arr[0][1] + middle_arr[0][1]) / 2) - (a_1[1] / (2 * a_2[1])),
        ]
        star_func_value = cur_func(star_x[0], star_x[1])
        print(a_1, a_2, star_x, star_func_value)
    return [star_x, star_func_value]


def gold_section(left_border_arr, least_func_value_arr, right_border_arr, gradient):

    print("\nЗОЛОТОЕ СЕЧЕНИЕ\n")
    alpha = right_border_arr[2] - left_border_arr[2]
    print(f"Gradient: {gradient}  \nAlpha: {alpha}\n")
    x1_lambda = left_border_arr[2] + 0.382 * alpha
    x2_lambda = left_border_arr[2] + 0.618 * alpha

    print(
        f"Lambdas : \nLeft border: {left_border_arr[2]}  \nX1_lambda: {x1_lambda} \
        \nMiddle: {least_func_value_arr[2]} \nX2_lambda: {x2_lambda} \
        \nRight_border: {right_border_arr[2]}\n"
    )

    x1_dot = [
        left_border_arr[0][0] - x1_lambda * gradient[0],
        left_border_arr[0][1] - x1_lambda * gradient[1],
    ]
    x2_dot = [
        left_border_arr[0][0] - x2_lambda * gradient[0],
        left_border_arr[0][1] - x2_lambda * gradient[1],
    ]

    x1_dot_func = cur_func(x1_dot[0], x1_dot[1])
    x2_dot_func = cur_func(x2_dot[0], x2_dot[1])

    print(
        f"Points: \nLeft: {left_border_arr[0], left_border_arr[1]} \nX1: {x1_dot, x1_dot_func} \
         \nMiddle: {least_func_value_arr[0], least_func_value_arr[1]} \nX2 : {x2_dot,x2_dot_func} \
         \nRight: {right_border_arr[0],right_border_arr[1]}"
    )

    while left_border_arr[1] > x1_dot_func and right_border_arr[1] > x2_dot_func:
        if left_border_arr[1] < right_border_arr[1]:
            right_border_arr = [x2_dot, x2_dot_func, x2_lambda]

        else:
            left_border_arr = [x1_dot, x1_dot_func, x1_lambda]

        alpha = right_border_arr[2] - left_border_arr[2]

        print(f"Gradient: {gradient}  \nAlpha: {alpha}\n")
        global tochnost_odn
        if alpha <= tochnost_odn:
            break

        x1_lambda = left_border_arr[2] + 0.382 * alpha
        x2_lambda = left_border_arr[2] + 0.618 * alpha

        print(
            f"Lambdas : \n Left border: {left_border_arr[2]}  \n X1_lambda: {x1_lambda} \
            \n Middle: {(right_border_arr[2] + left_border_arr[2])/2} \nX2_lambda: {x2_lambda} \
            \nRight_border: {right_border_arr[2]}\n"
        )

        x1_dot = [
            left_border_arr[0][0] - x1_lambda * gradient[0],
            left_border_arr[0][1] - x1_lambda * gradient[1],
        ]
        x2_dot = [
            left_border_arr[0][0] - x2_lambda * gradient[0],
            left_border_arr[0][1] - x2_lambda * gradient[1],
        ]

        x1_dot_func = cur_func(x1_dot[0], x1_dot[1])
        x2_dot_func = cur_func(x2_dot[0], x2_dot[1])

        print(
            f"Points: \n Left: {left_border_arr[0], left_border_arr[1]} \nX1: {x1_dot, x1_dot_func} \
             \nMiddle: {least_func_value_arr[0], least_func_value_arr[1]} \nX2 : {x2_dot,x2_dot_func} \
             \nRight: {right_border_arr[0],right_border_arr[1]}"
        )
    mindot = left_border_arr
    if mindot[1] > x1_dot_func:
        mindot = [x1_dot, x1_dot_func]
    if mindot[1] > x2_dot_func:
        mindot = [x2_dot, x2_dot_func]
    if mindot[1] > right_border_arr[1]:
        mindot = right_border_arr

    return mindot


def Pearson3(
    start_point, start_func_value, start_gradient, start_point_norm, start_norm_grad
):
    global criteriy_okonchania
    cur_point = start_point
    cur_func_value = start_func_value
    cur_gradient = start_gradient
    cur_point_norm = start_point_norm
    cur_norm_grad = start_norm_grad
    a = [1, 0, 0, 1]
    end = [100, 100]
    #cur_norm_grad > criteriy_okonchania
    while end[0] > criteriy_okonchania or end[1] > criteriy_okonchania:
        print(f"\n\nНОВАЯ ИТЕРАЦИЯ ДФП\n\n")
        dirrection = [
            a[0] * cur_gradient[0] + a[1] * cur_gradient[0],
            a[2] * cur_gradient[1] + a[3] * cur_gradient[1],
        ]
        print(f"Direction: {[-dirrection[0], -dirrection[1]]}\n")
        interval = Sven_algorithm(
            find_delta_lambda(cur_point_norm, find_norm(dirrection)),
            cur_point,
            cur_func_value,
            dirrection,
        )
        print(f"Sven result: {interval}\n")

        # (DSK_Pauela или gold_section)
        final_result = gold_section(interval[0], interval[1], interval[2], cur_gradient)
        print(f"Final: {final_result}")
        new_point = final_result[0]
        new_func_value = final_result[1]

        end = System_End_Criteria(new_point, new_func_value, cur_point, cur_func_value)
        print(f"Criteria: {end}\n")
        #cur_norm_grad < criteriy_okonchania
        if end[0] < criteriy_okonchania and end[1] < criteriy_okonchania:
            break
        new_grad = find_gradient(new_point, new_func_value)

        delta_x = [new_point[0] - cur_point[0], new_point[1] - cur_point[1]]
        delta_g = [new_grad[0] - cur_gradient[0], new_grad[1] - cur_gradient[1]]

        first_up = [
            delta_x[0] * delta_x[0],
            delta_x[0] * delta_x[1],
            delta_x[1] * delta_x[0],
            delta_x[1] * delta_x[1],
        ]
        first_down = delta_x[0] * delta_g[0] + delta_x[1] * delta_g[1]
        first_go = [
            first_up[0] / first_down,
            first_up[1] / first_down,
            first_up[2] / first_down,
            first_up[3] / first_down,
        ]

        second_up_one = [
            a[0] * delta_g[0] + a[1] * delta_g[1],
            a[2] * delta_g[0] + a[3] * delta_g[1],
        ]
        second_up_two = [
            second_up_one[0] * delta_g[0],
            second_up_one[0] * delta_g[1],
            second_up_one[1] * delta_g[0],
            second_up_one[1] * delta_g[1],
        ]
        second_up_three = [
            second_up_two[0] * a[0] + second_up_two[1] * a[2],
            second_up_two[0] * a[1] + second_up_two[1] * a[3],
            second_up_two[2] * a[0] + second_up_two[3] * a[2],
            second_up_two[2] * a[1] + second_up_two[3] * a[3],
        ]

        second_down_one = [
            delta_g[0] * a[0] + delta_g[1] * a[2],
            delta_g[0] * a[1] + delta_g[1] * a[3],
        ]
        second_down_two = (
            second_down_one[0] * delta_g[0] + second_down_one[1] * delta_g[1]
        )

        second_go = [
            second_up_three[0] / second_down_two,
            second_up_three[1] / second_down_two,
            second_up_three[2] / second_down_two,
            second_up_three[3] / second_down_two,
        ]

        a = [
            a[0] + first_go[0] - second_go[0],
            a[1] + first_go[1] - second_go[1],
            a[2] + first_go[2] - second_go[2],
            a[3] + first_go[3] - second_go[3],
        ]
        print(f"A: {a[0], a[1]}\n   {a[2], a[3]}\n")
        cur_point = final_result[0]
        cur_point_norm = find_norm(cur_point)
        cur_func_value = final_result[1]
        cur_gradient = find_gradient(cur_point, cur_func_value)
        cur_norm_grad = find_norm(cur_gradient)


if __name__ == "__main__":

    func_usage = 0

    Proizvodnie = 0.01
    parametr_Svena = 0.0001
    criteriy_okonchania = 0.001
    tochnost_odn = 0.0001
    h = Proizvodnie

    cur_point = [-1.2, 0]
    cur_norm_point = find_norm(cur_point)
    cur_func_value = cur_func(cur_point[0], cur_point[1])
    cur_gradient = find_gradient(cur_point, cur_func_value)
    cur_norm_gradient = find_norm(cur_gradient)

    Pearson3(cur_point, cur_func_value, cur_gradient, cur_norm_point, cur_norm_gradient)
    print(f"Колличество вызова функции: {func_usage}")

