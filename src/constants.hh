/***************************************************************************
                          constants.h  -  description
                             -------------------
    begin                : Thu Apr 11 2002
    copyright            : (C) 2002 by Mirela Andronescu
    email                : andrones@cs.ubc.ca
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef CONSTANTS_H
#define CONSTANTS_H

#define NONE            'N'         // no structure
#define HAIRP           'H'         // closes a hairpin loop
#define INTER           'I'         // closes an internal loop
#define MULTI           'M'         // closes a regular multi-loop

#define M_WM            'B'         // closes a regular partial multi-loop
#define M_WMv            'v'         // closes a regular partial multi-loop
#define M_WMp            'p'         // closes a regular partial multi-loop

#define FREE            'W'         // this base is free to be paired or unpaired to other base
#define LOOP            'V'         // closes a loop


#define P_P				'P'
#define P_PK			'k'
#define P_PL			'l'
#define P_PR			'r'
#define P_PM			'm'
#define P_PO			'o'
#define P_PfromL		'f'
#define P_PfromR		'g'
#define P_PfromM		'h'
#define P_PfromO		'i'
#define P_PLiloop		'j'
#define P_PLiloop5	'b'
#define P_PLmloop		'c'
#define P_PLmloop10	'e'
#define P_PLmloop01	'n'
#define P_PLmloop00	'a'
#define P_PRiloop		'q'
#define P_PRiloop5	's'
#define P_PRmloop		't'
#define P_PRmloop10	'u'
#define P_PRmloop01	'&'
#define P_PRmloop00	'9'
#define P_PMiloop		'w'
#define P_PMiloop5	'x'
#define P_PMmloop		'y'
#define P_PMmloop10	'0'
#define P_PMmloop01	'1'
#define P_PMmloop00	'8'
#define P_POiloop		'z'
#define P_POiloop5	'5'
#define P_POmloop		'+'
#define P_POmloop10	'-'
#define P_POmloop01	'='
#define P_POmloop00	'_'
#define P_WB			'*'
#define P_WBP			'^'
#define P_WP			'#'
#define P_WPP			'@'

#endif

