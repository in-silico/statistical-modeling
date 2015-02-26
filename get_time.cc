using namespace std;
#include <bits/stdc++.h>

typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::milliseconds milliseconds;

int main() {

  string line;
  while (getline(cin, line)) {
    if (line.size() == 0) continue;
    int exp = atoi(line.c_str());
    if (exp == -1)
      break;

    getline(cin, line);
    Clock::time_point t0 = Clock::now();
    getline(cin, line);
    Clock::time_point t1 = Clock::now();
    milliseconds ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
    std::cout << exp << " " <<  ms.count() << endl;
  }
  return 0;
}
