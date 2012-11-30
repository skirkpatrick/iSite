iSite: isite.cpp
	g++ -g isite.cpp -DNDEBUG -o iSite -Wall -pedantic -O3

debug: isite.cpp
	g++ -g isite.cpp -DDEBUG -o iSite -Wall -pedantic -O3

clean:
	rm -f iSite
