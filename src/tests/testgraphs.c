/* Unit tests for the graphs module.
 * Copyright (C) 2016-2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include "graphs.h"

void test2()
{
  matrix *dmx;
  NewMatrix(&dmx, 5, 5);
  G *graph;

  dmx->data[0][1] = dmx->data[1][0] = 0.3;
  dmx->data[0][2] = dmx->data[2][0] = 0.4;
  dmx->data[0][3] = dmx->data[3][0] = 0.7;
  dmx->data[0][4] = dmx->data[4][0] = 0.5;
  dmx->data[1][2] = dmx->data[2][1] = 0.9;
  dmx->data[1][3] = dmx->data[3][1] = 0.2;
  dmx->data[1][4] = dmx->data[4][1] = 0.2;
  dmx->data[2][3] = dmx->data[3][2] = 0.1;
  dmx->data[2][4] = dmx->data[4][2] = 0.6;
  dmx->data[3][4] = dmx->data[4][3] = 0.8;
  PrintMatrix(dmx);
  NewGraph(&graph);
  GenerateAdjMX(dmx, 0.8, &graph);
  PrintGraph(graph);
  DelMatrix(&dmx);
  DelGraph(&graph);
}


void test1()
{
  matrix *dmx;
  NewMatrix(&dmx, 4, 4);
  G *graph;

  dmx->data[0][1] = dmx->data[1][0] = 0.3;
  dmx->data[0][2] = dmx->data[2][0] = 0.4;
  dmx->data[0][3] = dmx->data[3][0] = 0.7;
  dmx->data[1][2] = dmx->data[2][1] = 0.9;
  dmx->data[1][3] = dmx->data[3][1] = 0.2;
  dmx->data[2][3] = dmx->data[3][2] = 0.1;
  PrintMatrix(dmx);
  NewGraph(&graph);
  GenerateAdjMX(dmx, 0.8, &graph);
  PrintGraph(graph);
  DelMatrix(&dmx);
  DelGraph(&graph);
}

int main(void)
{
  test1();
}
