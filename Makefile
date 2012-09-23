iSite: isite.cpp
	g++ -g isite.cpp -DNDEBUG -o iSite

debug: isite.cpp
	g++ -g isite.cpp -DDEBUG -o iSite

clean:
	rm -f iSite
