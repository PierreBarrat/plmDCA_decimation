per compilare ho fatto:
g++ do_montecarlo.cpp graph1.cpp graph2.cpp graph3.cpp graphs.cpp gasdev.cpp -o do_montecarlo.out

per eseguire:
./do_montecarlo.out -n 53 -q 21 -m 100 < example_pars.txt



!!!!!!!!!!!!!!!attenzione
ho cambiato la riga 137 e 138 di graph2.cpp
1e-3 -> 1e^3 e' okok


anche qui:
do_montecarlo.out: graph2.cpp:155: void Graph2::read(std::istream&): Assertion `fabs(sumh) < 1e-3' failed.

