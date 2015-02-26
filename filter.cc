using namespace std;
#include <bits/stdc++.h>

int main() {
  vector<int> seen(41, -1);
  string line;
  while (getline(cin, line)) {
    stringstream ss(line);
    int exp, time;
    ss >> exp >> time;
    if (seen[exp] >= 0)
      continue;
    seen[exp] = time;
  }

  for (int i = 1; i < seen.size(); ++i) {
    if (seen[i] < 0)
      cerr << "missing " << i << endl;
    else
      cout << i << " " << seen[i] << endl;
  }

  return 0;
}
