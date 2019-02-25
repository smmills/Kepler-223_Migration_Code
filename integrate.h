//  Written by Cezary Migazewski - Last Edit 1/1/16

void integrate ( my_float Y1[], my_float Y2[], my_float time1, my_float time2, struct params_str params )
{

// rk4(&equations, Y1, Y2, time1, time2-time1, params);
// rk2(&equations, Y1, Y2, time1, time2-time1, params);
// lobattoIIIA_2nd(&equations, Y1, Y2, time1, time2-time1, params);
// lobattoIIIA_4th(&equations, Y1, Y2, time1, time2-time1, params);
// lobattoIIIAs4(&equations, Y1, Y2, time1, time2-time1, params);
 lobattoIIIAs5(&equations, Y1, Y2, time1, time2-time1, params);
// lobattoIIICs5(&equations, Y1, Y2, time1, time2-time1, params);

}
