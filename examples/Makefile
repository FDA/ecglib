all:	getdbannotations twavedelineator 

getdbannotations:
	g++ -std=c++11 -o getdbannotations getdbannotations.cpp -lboost_program_options -lboost_filesystem -lboost_system -lwfdb
twavedelineator:
	g++ -std=c++11 -o twavedelineatorphysionet twavedelineatorphysionet.cpp filters/butterworth/butter.cpp -lecglib-core -lecglib-delineators-twave -lboost_program_options -lboost_filesystem -lboost_system -lwfdb
clean:
	rm getdbannotations twavedelineatorphysionet
