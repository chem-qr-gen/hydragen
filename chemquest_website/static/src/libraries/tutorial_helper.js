import Shepherd from "shepherd.js";

function initiateTutorial() {
  setTimeout(applyFilter(null), 0);
  welcomeText();

  $("#overlay").on("click", async () => {
    $("#overlay").on("click", async () => {
    })
    $("#overlay").toggleClass("show");
    startTour();
  });
}

function applyFilter(element) {
  var clipPath = `polygon(
     0% 0%,
     0% 100%,
     100% 100%,
     100% 0%
  )`;
  if (element !== null) {
    // const mainRect = $("#main-container")[0].getBoundingClientRect();
    const rect = element[0].getBoundingClientRect();
    // console.log(mainRect.top, mainRect.right, mainRect.bottom, mainRect.left);
    // console.log(rect.top, rect.right, rect.bottom, rect.left);

    const padding = 10;

    const top = rect.top - padding;
    const bottom = rect.bottom + padding;
    const left = rect.left - padding;
    const right = rect.right + padding;

    clipPath = `polygon(
      0% 0%,
      0% 100%,
      ${left}px 100%,
      ${left}px ${top}px,
      ${right}px ${top}px,
      ${right}px ${bottom}px,
      ${left}px ${bottom}px,
      ${left}px 100%,
      100% 100%,
      100% 0%,
      0% 0%
    )`;
  }

  $("#overlay").css("clip-path", clipPath);
  $("#overlay").css("-webkit-clip-path", clipPath);

}

function welcomeText() {
  moveTB1("40%", "50%");
  moveTB2("50%", "50%")
  $("#tutorial-text2").addClass('flash');
}

function moveTB1(top, left) {
  if (top !== null) {
    $("#tutorial-text-box1").css("top", top);
  }
  if (left !== null) {
    $("#tutorial-text-box1").css("left", left);
  }
}

function moveTB2(top, left) {
  if (top !== null) {
    $("#tutorial-text-box2").css("top", top);
  }
  if (left !== null) {
    $("#tutorial-text-box2").css("left", left);
  }

}

function startTour() {
  var tour = new Shepherd.Tour({
    useModalOverlay: true,
    defaultStepOptions: {
      classes: "custom-shepherd-step",
      cancelIcon: {
        enabled: true
      },
      scrollTo: true,
      // cancelIcon: false
    }
  });

  tour.addStep({
    text: "Welcome to Hydragen! Hygragen is an app where you can learn mass spectrometry in a fun and interactive way.",
    buttons: [
      {
        action() {
          return this.back();
        },
        classes: 'shepherd-button-secondary',
        text: 'Back'
      },
      {
        action() {
          return this.next();
        },
        classes: "custom-shepherd-button",
        text: 'Next'
      }
    ],
    id: ''
  });

  tour.addStep({
    text: "Here's the mass spectrometry graph.\nIf you are stuck, click on a bar in the chart for a hint!",
    attachTo: {
      element: $("#msChart")[0],
      on: 'right'
    },
    buttons: [
      {
        action() {
          return this.back();
        },
        classes: 'shepherd-button-secondary',
        text: 'Back'
      },
      {
        action() {
          return this.next();
        },
        classes: "custom-shepherd-button",
        text: 'Next'
      }
    ],
    id: ''
  });

  tour.addStep({
    text: "Here shows how many hints you have left. You have 3 hints for every question.",
    attachTo: {
      element: $("#hint-display")[0],
      on: 'bottom-start'
    },
    buttons: [
      {
        action() {
          return this.back();
        },
        classes: 'shepherd-button-secondary',
        text: 'Back'
      },
      {
        action() {
          return this.next();
        },
        classes: "custom-shepherd-button",
        text: 'Next'
      }
    ],
    id: ''
  });

  tour.addStep({
    text: "Select the structure you think would produce the MS chart.\nYou have multiple attempts, but you'll get a lower score with each attempt.",
    attachTo: {
      element: $("#answer")[0],
      on: 'left'
    },
    buttons: [
      {
        action() {
          return this.back();
        },
        classes: 'shepherd-button-secondary',
        text: 'Back'
      },
      {
        action() {
          return this.next();
        },
        classes: "custom-shepherd-button",
        text: 'Next'
      }
    ],
    id: ''
  });

  tour.addStep({
    text: "Click on option to select it.",
    attachTo: {
      element: $("#radio-opt0")[0],
      on: 'left'
    },
    buttons: [
      {
        action() {
          return this.back();
        },
        classes: 'shepherd-button-secondary',
        text: 'Back'
      },
      {
        action() {
          return this.next();
        },
        classes: "custom-shepherd-button",
        text: 'Next'
      }
    ],
    id: ''
  });

  tour.addStep({
    text: 'Click here to view the periodic table.',
    attachTo: {
      element: $("#switch-graph")[0],
      on: 'top'
    },
    buttons: [
      {
        action() {
          return this.back();
        },
        classes: 'shepherd-button-secondary',
        text: 'Back'
      },
      {
        action() {
          return this.next();
        },
        classes: "custom-shepherd-button",
        text: 'Next'
      }
    ],
    id: ''
  });

  tour.addStep({
    text: 'Click here to submit the answer you have chosen. Think carefully before you submit!',
    attachTo: {
      element: $("#submit")[0],
      on: 'top'
    },
    buttons: [
      {
        action() {
          return this.back();
        },
        classes: 'shepherd-button-secondary',
        text: 'Back'
      },
      {
        action() {
          return this.next();
        },
        classes: "custom-shepherd-button",
        text: 'Next'
      }
    ],
    id: ''
  });

  tour.addStep({
    text: 'If you are stuck, click here to skip the question.',
    attachTo: {
      element: $("#next")[0],
      on: 'top'
    },
    buttons: [
      {
        action() {
          return this.back();
        },
        classes: 'shepherd-button-secondary',
        text: 'Back'
      },
      {
        action() {
          return this.next();
        },
        classes: "custom-shepherd-button",
        text: 'Next'
      }
    ],
    id: ''
  });

  tour.addStep({
    text: "Good Luck Have Fun",
    buttons: [
      {
        action() {
          return this.back();
        },
        classes: 'shepherd-button-secondary',
        text: 'Back'
      },
      {
        action() {
          return this.next();
        },
        classes: "custom-shepherd-button",
        text: "Let\'s GO!"
      }
    ],
    id: ''
  });

  tour.start();
}

export {
  initiateTutorial
}
