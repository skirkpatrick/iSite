iSite: isite.cpp
	g++ -g isite.cpp -o iSite

debug: isite.cpp
	g++ -g isite.cpp -DDEBUG -o iSite

clean:
	rm -f iSite
